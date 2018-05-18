use strict;

require Utils;

my (%tasks, $goMeta, $goClosure);
my $goSrc   = "ftp://ftp.geneontology.org/go/ontology/go.obo";
my $defEC   = "UNK,ND,P,IC,IEA,NAS,IRD,IKR,IBD,IBA,RCA,IGC,ISS,ISA,ISM,ISO,TAS,EXP,IEP,IPI,IMP,IGI,IDA";

## GO Relation types
## http://www.geneontology.org/page/ontology-relations
my $defRel  = "regulates,positively_regulates,negatively_regulates,intersection_of,part_of,is_a";

our ($bar, $args);

$args->{rellevels} ||= $defRel;
$args->{eclevels}  ||= $defEC;

     
## Stash deduplicated file store - not on all systems
my $stashMeta = {
    Authority  => "GeneOntology",
    FileType   => "AnnotatedMatrix",
    Format     => "MatrixMarket",
};

my $goUrl = {
    GeneOntology   => 'http://amigo.geneontology.org/amigo/term/%s',
};

my @stndGoMeta = qw(Name Namespace Description Parents Synonyms);
my %defGoDef = ( Name => "A short name associated with the GO accession, eg 'actin binding'",
    Namespace => "The high-level GO tree the term belongs to",
    Description => "A longer definition of a GO term, taken from the 'def' field",
    Parents => "The parent terms associated with the GO term, either through is_a or part_of, concatenated with '|'",
    Synonyms => "Alternate (unofficial) short names for the GO term, concatenated with '|'");

sub make_go_matrix {
    ## Make a matrix between objects (eg genes) and GO terms
    my ($links,   # Hash keyed to object, values array of [GO,EC,GO,EC,...]
        $objMeta, # Metadata for the object IDs
        $trg,     # The matrix file path to be created
        $fbits,   # File name components
        $fmeta,   # File metadata
        $nsUrl,   # URL templates as hash, keyed to namespace
        $authL,   # "Long" authority name, for matrix metadata
        $citeCB,  # Function reference to generate MTX citation block
        $taxDat,  # Taxonomy information built via Utils.pm
        $objMcol, # Metadata columns to report for the objects
        $objMdef, # Descriptions for the object metadata columns
        ) = @_;

    my ($lvls, $scH) = &_factorize_levels('eclevels', 'upper');
    &go_meta();
    ## Make a full child-to-ancestor lookup structure:
    my $closure = &go_transitive_closure(['is_a','part_of']);
    my %expand;
    while (my ($kid, $pars) = each %{$closure}) {
        $expand{$kid} = [ $kid, keys %{$pars} ];
    }


    my (%objAssignments, %ontOrder); ## Was %genes and %ontMeta
    my ($nznum, $nont) = (0,0);
    warn "     ... structuring GO links through transitive closure ...\n";
    my @objIds = sort keys %{$links};
    foreach my $oId (@objIds) {
        my @goEc = @{$links->{$oId}};
        my $targ = $objAssignments{$oId} = {};
        for (my $i = 0; $i < $#goEc; $i += 2) {
            my ($kid, $ec) = ($goEc[$i], $goEc[$i+1]);
            my @allGo = @{$expand{$kid} || []};
            if ($#allGo == -1) {
                &msg("[DataError] - Could not find ancestors for $kid")
                    unless ($tasks{"CarpNoParents-$kid"}++);
                next;
            }
            my $sc = $scH->{ uc($ec) };
            unless ($sc) {
                unless (defined $sc) {
                    &err("Unknown Evidence code '$ec'",
                         "This code was not in -eclevels, it will be ignored");
                    $scH->{ uc($ec) } = 0;
                }
                next;
            }
            
            foreach my $gid (@allGo) {
                my $gm = $ontOrder{$gid} ||= {
                    ## Just to track order of GO column indices
                    order => ++$nont,
                };
                my $gnum = $gm->{order};
                # Capture orderID/score pairs
                if (!$targ->{$gnum}) {
                    ## First time we've seen this GO term for this object
                    $nznum++;
                    ## Just set the score as-is
                    $targ->{$gnum} = $sc;
                } elsif ($targ->{$gnum} < $sc) {
                    ## We have seen this GO term before, but with a
                    ## "worse" evidence code. "Upgrade" the EC to
                    ## reflect this better evidence.
                    $targ->{$gnum} = $sc;
                    ## This can/will happen becuase of inheritance
                }
            }
        }
    }
    
    warn "     ... writing matrix market file ...\n";
    my @ontIds = sort { $ontOrder{$a}{order} 
                        <=> $ontOrder{$b}{order} } keys %ontOrder;
    my ($rnum, $cnum) = ($#objIds + 1, $#ontIds + 1);
    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write GO ontology", $tmp, $!);

    my ($auth, $mod, $nsi, $nsj, $type) = 
        map { $fbits->{$_} } qw(auth mod ns1 ns2 type);

    print MTX &_initial_mtx_block
        ("Ontology", $rnum, $cnum, $nznum, "$auth $mod $nsi-to-$nsj Ontology",
         "$mod $nsi $nsj ontology, as assigned by $auth",
         "Scores are factors representing the evidence code for the assignment",
         $nsi, $nsj);

    my %nsu =( %{$goUrl}, %{$nsUrl} );

    print MTX &_dim_block({
        %{$fmeta},
        RowDim    => $nsi,
        RowUrl    => $nsu{$nsi},
        ColDim    => $nsj,
        ColUrl    => $nsu{$nsj}, 
        Authority => $authL });

    print MTX &{$citeCB}() if ($citeCB);
    print MTX &_species_MTX( $taxDat->{'scientific name'} ) if ($taxDat);

    print MTX &_filter_block({
        AutoFilterComment => "NOTE: The standard automatic filters for this matrix are designed with enrichment analysis in mind. If you are using this matrix to recover GO terms for genes (or vice versa), you should \$reset() the matrix, and optionally re-apply only the evidence code based filters.",
        MinColCount => "7 ## Terms with few genes assigned to them struggle to reach statistical significance",
        MaxColCount => "10% ## Terms covering a large fraction of the genome are rarely informative. Note this will exclude many parent terms",
        MinRowCount => "2 ## Genes with few terms assigned to them will not bring much insight to the analysis",
        TossLevel => "[,][ND,IC,P] ## ND=No Data, IC=Inferred by Curator, P is an ancient evidence code"
                             });

    my @counts;
    foreach my $h (values %objAssignments) {
        my @u   = values %{$h};
        map { $counts[$_-1]++ } @u;
    }
    print MTX &_factor_map_block
        ( $lvls, $scH, "Evidence Codes",
          ["Lower factor values generally represent lower confidence evidence",
           "http://geneontology.org/page/guide-go-evidence-codes"],
          {map { $lvls->[$_] => $counts[$_] || 0 } (0..$#{$lvls}) });

    print MTX &_rowcol_meta_comment_block( {%{$objMdef}, %defGoDef, Description => "Descriptive text (either for Entrez ID or GO term)"});

    ## Row metadata
    print MTX &_generic_meta_block(\@objIds, 'Row', $objMeta, $objMcol);
    print MTX "% $bar\n";

    ## Column metadata
    print MTX &_generic_meta_block(\@ontIds, 'Col', $goMeta, \@stndGoMeta);

    print MTX &_triple_header_block( $rnum, $cnum, $nznum );
    for my $i (0..$#objIds) {
        while (my ($j, $sc) = each %{$objAssignments{$objIds[$i]}}) {
            printf(MTX "%d %d %d\n", $i+1, $j, $sc);
        }
    }
    close MTX;
    rename($tmp, $trg);
    &msg("Generated $nsj $type", $trg);
    &post_process( %{$fbits}, meta => $fmeta );
    return $trg;
}

sub make_go_metadata_file {
    ## Simple TSV file of GeneOntology metadata

    my $src   = $goSrc;
    my $file  = &fetch_url($src);
    
    ## Get the version
    my $halt  = 0;
    open(GOFILE, "<$file") || 
        &death("Failed to read metadata file", $file, $!);
    my $fVers;
    while (<GOFILE>) {
        &death("Did not find version line!", $file) if (++$halt > 100);
        if (/data-version: releases\/(\d{4}-\d{2}-\d{2})/) {
            $fVers = $1;
            last;
        }
    }
    close GOFILE;

    my $type  = "Metadata";
    my $cont  = "GeneOntology"; ## Both the content and authority
    my %fbits = (type => $type, sfx => 'tsv',
                 ns1  => $cont,
                 auth => $cont,  vers => $fVers);
    my $meta = {
        MatrixType => $fbits{type},
        Modifier   => $fbits{mod},
        Species    => undef,
        Authority  => $cont,
        FileType   => "TSV",
        Version    => $fVers,
        'Format'   => "TSV",
        Source     => &ftp_url($src),
    };
    my $fmeta = { %{$stashMeta}, %{$meta} };
    my $trg   = &primary_path(%fbits);

    ## ignoring has_part : inverse of part_of

    unless (&output_needs_creation($trg)) {
        unless ($tasks{GoMeta}++) {
            &msg("Using existing Metadata file:", $trg);
            &post_process( %fbits, meta => $fmeta );
        }
        return $trg;
    }

    &msg("Parsing basic GeneOntology metadata");
    my $tmp = "$trg.tmp";
    open (METAF, ">$tmp") || &death("Failed to write metadata file", $tmp, $!);
    my @head = qw(id name Parents is_obsolete namespace is_a part_of intersection_of regulates positively_regulates negatively_regulates synonym def);
    my $hMap = {
        id          => "ID",
        name        => "Name",
        namespace   => "Namespace",
        def         => "Description",
        synonym     => "Synonyms",
        is_obsolete => "Obsolete",
    };
    print METAF join("\t", map { $hMap->{$_} || $_ } @head)."\n";

    open(GOFILE, "<$file") || &death("Failed to read metadata file", $file, $!);
    my $rec = { num => 0, now => "HEAD", dat => {} };
    while (<GOFILE>) {
        s/[\n\r]+$//;
        if (/^\[(\S+)\]/) {
            &_go_meta_record( $rec, *METAF, \@head, $1 );
        } elsif (/^(\S+):\s+(.+?)\s*$/) {
            my ($k, $v) = ($1, $2);
            if ($v =~ /^\"(.+)\"/) {
                $v = $1;
                $v =~ s/\\"/"/g;
            } else {
                $v =~ s/\s+!.+//; # Remove ! comments
            }
            if ($k eq 'relationship') {
                if ($v =~ /^(\S+)\s+(GO:\d{7})$/) {
                    ($k, $v) = ($1, $2);
                } else {
                    &err("Weird GO line: $_");
                    next;
                }
            } elsif ($k eq 'intersection_of') {
                if ($v =~ / (GO:\d{7})$/) { $v = $1 }
            }
            push @{$rec->{dat}{$k}}, $v;
        }
    }
    &_go_meta_record( $rec, *METAF, \@head, "end" );
    close GOFILE;
    close METAF;
    rename($tmp, $trg);
    &msg("Generated GeneOntology metadata file", $trg);
    $tasks{GoMeta}++;
    &post_process( %fbits, meta => $fmeta );
    return $trg;
}

sub _go_meta_record {
    my ($rec, $fh, $head, $nextK) = @_;
    my $k = $rec->{now};
    my $d = $rec->{dat};
    if ($k eq 'Term') {
        if (my $id = $d->{id}) {
            my @row;
            ## Build the parent array from is_a and part_of
            my %u;
            foreach my $t ('is_a', 'part_of') {
                if (my $p = $d->{$t}) {
                    map { $u{$_} = 1 } @{$p};
                }
            }
            delete $u{""};
            my @par = sort keys %u;
            $d->{Parents} = \@par unless ($#par == -1);
            foreach my $c (@{$head}) {
                if (my $vs = $d->{$c}) {
                    map { s/\|/_/g } @{$vs};
                    my %u = map { $_ => 1 } @{$vs};
                    $vs = [ sort keys %u ];
                    push @row, join('|', @{$vs});
                } else {
                    push @row, "";
                }
            }
            print $fh join("\t", @row)."\n";
        } else {
            warn Dumper($d);
            &death("Failed to parse GO ID from record")
        }
    }
    $rec->{dat} = {};
    $rec->{now} = $nextK;
}

sub go_meta {
    return $goMeta if ($goMeta);
    $goMeta = {};
    my $mdFile = &make_go_metadata_file();
    my ($lvls, $scH) = &_factorize_levels('rellevels');
    my $nLvl = $#{$lvls};
    open(METAF, "<$mdFile") || &death("Failed to read metadata TSV",
                                      $mdFile, $!);
    my $head = <METAF>;
    $head =~ s/[\n\r]+$//;
    my @cols = split(/\t/, $head);
    while (<METAF>) {
        s/[\n\r]+$//;
        my @row  = split(/\t/);
        my %data = map { $cols[$_] => $row[$_] } (0..$#cols);
        my $gid  = $data{ID};

        if ($goMeta->{$gid}) {
            &err("Multiple entries for GeneID == $gid");
        } else {
            $goMeta->{$gid} = \%data;
            ## Set up parentage hash
            for my $i (0..$#{$lvls}) {
                my $targ = $goMeta->{PARENTAGE}{$gid} ||= {};
                if (my $pars = $data{$lvls->[$i]}) {
                    my $sc   = $i + 1;
                    foreach my $par (split(/\|/, $pars)) {
                        $targ->{$par} = $sc
                            if (!$targ->{$par} || $targ->{$par} < $sc);
                    }
                }
            }
        }
    }
    close METAF;
    return $goMeta;
}

sub go_transitive_closure {
    my $lvlReq = shift;
    
    if ($lvlReq) {
        ## A request to limit the closure to only particular relationships:
        $lvlReq = split(/\s*,\s*/, $lvlReq) unless (ref($lvlReq));
        $lvlReq = { map { lc($_ || "") => 1 } @{$lvlReq} };
        delete $lvlReq->{""};
        &msg("Restricting transitive GO closure to: ".
             join(", ", sort keys %{$lvlReq}));
    } elsif ($goClosure) {
        return $goClosure;
    }

    &go_meta();
    my ($lvls, $scH) = &_factorize_levels('rellevels');
    ## http://www.geneontology.org/page/ontology-relations
    my $transitive = {
        "is_a" => {
            "is_a"      => "is_a",
            "part_of"   => "part_of",
            "regulates" => "regulates",
        },
        "part_of" => {
            "part_of"   => "part_of",
            "is_a"     => "part_of",
        },
        "regulates" => {
            "is_a"      => "regulates",
            "part_of"   => "regulates",
        },
        "positively_regulates" => {
            "is_a"      => "positively_regulates",
            "part_of"   => "regulates",
        },
        "negatively_regulates" => {
            "is_a"      => "negatively_regulates",
            "part_of"   => "regulates",
        },
    };
    ## Normalize as 2D score array:
    my @tranSC;
    while (my ($e1, $e2H) = each %{$transitive}) {
        next if ($lvlReq && !$lvlReq->{lc($e1)});
        my $sc1 = $scH->{uc($e1)};
        while (my ($e2, $e3) = each %{$e2H}) {
            next if ($lvlReq && !$lvlReq->{lc($e2)});
            $tranSC[$sc1][ $scH->{uc($e2)} ] = $scH->{uc($e3)};
        }
    }

    ## Begin by initializing the structure with the direct connections
    my $gc = {};
    foreach my $gid (keys %{$goMeta}) {
        next if ($goMeta->{$gid}{Obsolete});
        $gc->{$gid} = { %{$goMeta->{PARENTAGE}{$gid} || {}} };
        # die Dumper($gc);
    }

    my @nodes     = sort keys %{$gc};
    my $iteration = 1;
    my $perturbed = { map { $_ => 1 } @nodes };
    my $changed   = $#nodes + 1;
    &msg("Creating transitive closure for GO");
    while ($changed) {
        warn sprintf("  Iteration %d: %d changes\n", $iteration++, $changed);
        $changed = 0;
        my $ptb = {};
        foreach my $kid (@nodes) { # Kid ID
            foreach my $pid (keys %{$gc->{$kid}}) { # Parent ID
                # next unless ($perturbed->{$pid});
                my $sc1 = $gc->{$kid}{$pid};
                foreach my $gid (keys %{$gc->{$pid}}) {
                    my $sc2 = $gc->{$pid}{$gid};
                    # warn "$kid -[$sc1]-> $pid -[$sc2]-> $gid\n";
                    if (my $sc3 = $tranSC[$sc1][$sc2]) {
                        ## There is a transitive path through the parent
                        if (!$gc->{$kid}{$gid} || $gc->{$kid}{$gid} < $sc3) {
                            ## And that path is new, or better
                            $gc->{$kid}{$gid} = $sc3;
                            $changed++;
                            $ptb->{$kid}++;
                        }
                    }
                }
            }
        }
        $perturbed = $ptb;
    }
    ## Set the global value so long as we did not filter by level:
    $goClosure = $gc unless ($lvlReq);
    return $gc;
}

sub go_hierarchy {
    my $type  = shift || "Hierarchy";
    my $src   = $goSrc;
    my $fVers = &_datestamp_for_file(&make_go_metadata_file());
    my $cont  = "GeneOntology";    
    my ($nsi, $nsj) = ("Child", "Parent");
    my %fbits = (type => $type,  mod  => $cont,
                 ns1  => $nsi,   ns2  => $nsj,
                 auth => $cont,  vers => $fVers );

    my $meta = {
        MatrixType => $fbits{type},
        Modifier   => $fbits{mod},
        Source     => &ftp_url($src),
        Authority  => $cont,
        Namespace  => [$nsi, $nsj],
    };
    my $fmeta = { %{$stashMeta}, %{$meta} };
    my $trg   = &primary_path(%fbits);
    unless (&output_needs_creation($trg)) {
        unless ($tasks{"GO-$type"}++) {
            &msg("Keeping existing $cont $type file:", $trg);
            &post_process( %fbits, meta => $fmeta );
        }
        return $trg;
    }
    &go_meta();
    ## The Hierarchy will have only direct parent-child links. The
    ## Closure will include all transitively closed paths from a child
    ## to any valid ancestor.
    my $obj = $type eq 'Closure' ? 
        &go_transitive_closure() : $goMeta->{PARENTAGE};
    my ($lvls, $scH) = &_factorize_levels('rellevels');
    my @gids;
    foreach my $gid (sort keys %{ $obj }) {
        push @gids, $gid unless ($goMeta->{$gid}{Obsolete});
    }
    my $rnum  = $#gids + 1;
    my $nznum = 0;
    my %counts;
    foreach my $h (values %{ $obj }) {
        my @u   = values %{$h};
        $nznum += $#u + 1;
        map { $counts{ $lvls->[$_-1]}++ } @u;
    }
    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write $cont $type", $tmp, $!);
    print MTX &_initial_mtx_block
        ( $type, $rnum, $rnum, $nznum, "$cont $nsi-to-$nsj $type",
          "$nsi-$nsj relationships between $cont graph nodes",
          "Node relationship types, as factor levels",
          $nsi, $nsj);

    print MTX &_dim_block({
        %{$fmeta},
        RowDim     => $nsi,
        RowUrl     => $goUrl->{$cont},
        ColDim     => $nsj,
        ColUrl     => $goUrl->{$cont} });

    if ($type eq 'Closure') {
        print MTX &_mtx_comment_block("This matrix represents a full transitive closure of the GO hierarchy. Rows are child nodes, columns represent all direct or indirect ancestors (parents, grandparents, etc). Some ancestors might be accessible through multiple paths; If so, the relationship type will be the 'best' one (higher score value in the relationship factor). For direct relationships only, see the Hierarchy matrix.");
    } else {
        print MTX &_mtx_comment_block("This matrix represents direct $cont parentage relationships between the rows (children) and columns (parents). For full transitive parentage, see the Closure matrix.");
    }

    print MTX &_go_citation_MTX();
    print MTX &_filter_block({
        KeepLevel => "[,][part_of,is_a] ## 'Part Of' and 'Is A' relationships are the only ones that reliably allow inheritance graph building."
                             });

    print MTX &_factor_map_block
        ( $lvls, $scH, "Ontology Relations",
          ["The kinds of relationships that structure the GO graph",
           "Only is_a and part_of are reliable for transitive graph following",
           "http://www.geneontology.org/page/ontology-relations"],
        \%counts);

    print MTX &_rowcol_meta_comment_block(\%defGoDef);
    ## Row metadata
    print MTX &_generic_meta_block(\@gids, 'Row', $goMeta, \@stndGoMeta);

    ## Column metadata
    print MTX "% $bar\n";
    print MTX &_generic_meta_block(\@gids, 'Col', $goMeta, \@stndGoMeta);

    print MTX &_triple_header_block( $rnum, $rnum, $nznum );
    my $n2i = { map { $gids[$_] => $_+1 } (0..$#gids) };
    for my $i (0..$#gids) {
        while (my ($par, $sc) = each %{$obj->{$gids[$i]}}) {
            ## die "$gids[$i] + $par = $sc" unless ($n2i->{$par});
            printf(MTX "%d %d %d\n", $i+1, $n2i->{$par}, $sc);
        }
    }
    close MTX;
    rename($tmp, $trg);
    &msg("Generated $cont $nsi to $nsj $type", $trg);
    &post_process( %fbits, meta => $fmeta );
    return $trg;
}

sub _go_citation_MTX {
    my $fVers = &_datestamp_for_file(&make_go_metadata_file());
  
    return &_default_parameter( "Citation", "The Gene Ontology Consortium. Gene Ontology Consortium: going forward. (2015) Nucl Acids Res 43 Database issue D1049–D1056. Downloaded $fVers.")."%\n";
}
