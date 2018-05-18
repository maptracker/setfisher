#!/usr/bin/perl -w

use strict;
my $scriptDir;
our $defTmp = "/tmp/ensemblGeneSets";
my $defFtp  = "ftp.ensembl.org";

our $defaultArgs = {
    ftp       => $defFtp,
    dir       => ".",
};

my $revision = '1';
my $revNotes = {
    '0' => "# Beta code, still under development",
    '1' => "# Fix protein description",
};

## XREF DBs to take "as is"
my $asIsXrefs = {
    UniParc => 1,
    PDB => 1,
    MEROPS => 1,
    MetaCyc => 1,
    UniPathway => 1,
};

## XREF DBs we will quietly ignore
my $ignoreXrefs = {
    Clone_based_ensembl_transcript => 1,
    Clone_based_vega_transcript    => 1,
    EntrezGene_trans_name          => 1,
    'Uniprot/SPTREMBL'             => 1,
    miRBase_trans_name             => 1,
    EMBL            => 1,
    RGD_trans_name  => 1,
    RFAM_trans_name => 1,
    MGI_trans_name  => 1,
    RNAcentral      => 1,
    protein_id      => 1,
    Reactome        => 1,
    ChEMBL          => 1, # Anything of value here?
};

## Primary tags that don't seem to have useful info:
my $ignorePrimaryTag = {
    source          => 1,
    exon            => 1,
    misc_feature    => 1, # Appears to just be contig footprints?
    STS             => 1,
};

## Secondary tags that we are ignoring
my $ignoreFeatTags = {
    translation     => 1,
    codon_start     => 1,

};


BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}
use lib "$scriptDir";
require Utils;
require GoMetadata;

our ($args, $outDir, $clobber, $tmpDir, $maxAbst, $bar);
my (%tasks, $geneData, $ftpDir, $globalMeta);

my $taxDat     = &extract_taxa_info( $args->{species} || $args->{taxa} );
my $specID     = $args->{name} || &species_id_for_tax( $taxDat );
my $sciName    = $args->{name} || &species_scientific_name( $taxDat );
my $vers       = $args->{vers} || $args->{version};
if (!$vers) {
    $taxDat->{error} ||= "No -version provided";
} elsif ($vers =~ /^(release-)?(\d+)$/) {
    ## Allow it to be specified as 'release-92' as well as '92'
    $vers = $2;
} else {
    $taxDat->{error} ||= "-version is not an integer";
}

my $versToken = "release-$vers"; 
my $ftpSpec   = lc($sciName);    # homo_sapiens
$ftpSpec      =~ s/ /_/g;
my $auth      = "Ensembl";
my $authLong  = "$auth ## Major genome annotation effort by the EBI - https://www.ensembl.org/info/about/";

## Stash deduplicated file store - not on all systems
my $stashMeta = {
    Authority  => $auth,
    'Version'  => $vers,
    MatrixType => "Map",
    FileType   => "AnnotatedMatrix",
    Format     => "MatrixMarket",
};

## Identify and download EMBL files:
my $manifest = &fetch_all_files();

if ($taxDat->{error} || $args->{h} || $args->{help}) {
    warn "
Usage:

This program will generate gene sets for the SetFisher enrichment
package. It will recover information from the Ensembl FTP site based on
a species identifier and version number you provide.

Required Argument:

  -species The species name you wish to extract. This should match the
           Ensembl species name.

  -version The Ensembl version number. See here:

           ftp://ftp.ensembl.org/pub/

           ... available versions will start with 'release-'

Optional Arguments:

      -ftp Default '$defFtp'. URL to Ensembl's FTP site

   -tmpdir Default '$defTmp'. Directory holding downloaded files

      -dir Default '.'.  Output directory that will contain generated
           matrix and metadata files

  -clobber Default 0, which will preserve any already-generated
           files. Set to 1 to regenerate matrix files, and any value
           greater than 1 to also re-download base files from Entrez.

     -name A name to use for the files. By default this will be the
           Ensembl species name.


";
    warn "\n[!!] Failed to fulfill your request:\n     $taxDat->{error}\n\n"
        if ($taxDat->{error});
    exit;
}

my $ensSearch  = 'http://www.ensembl.org/Multi/Search/Results?q=%s';
my $nsUrl = {
    Symbol         => 'https://www.ncbi.nlm.nih.gov/gene/?term=%s%%5Bsym%%5D',
    EntrezGene     => 'https://www.ncbi.nlm.nih.gov/gene/%s', # Integer IDs
    LocusLink      => 'https://www.ncbi.nlm.nih.gov/gene/?term=%s', # For LOC###
    PubMed         => 'https://www.ncbi.nlm.nih.gov/pubmed/%s',
    GeneOntology   => 'http://amigo.geneontology.org/amigo/term/%s',
    RefSeqRNA      => 'https://www.ncbi.nlm.nih.gov/nuccore/%s',
    RefSeqProtein  => 'https://www.ncbi.nlm.nih.gov/protein/%s',
    Orthologue     => 'https://www.ncbi.nlm.nih.gov/taxonomy/?term=%s',
    EnsemblGene    => $ensSearch,
    EnsemblRNA     => $ensSearch,
    EnsemblProtein => $ensSearch,
};
my %defColDef = ( Symbol => "Official Ensembl Gene symbol, if available, otherwise may be a sequence accession",
                  Description => "Short descriptive text for an Ensembl accession" );

&go_hierarchy("Hierarchy");
&go_hierarchy("Closure");
&go_ontology();

sub basic_links {
    ## Gene <-> RNA, RNA <-> Protein
    my %fbitCommon = (type => "Map",  mod  => $specID,
                      ## ns1  => $nsi,   ns2  => $nsj,
                      auth => $auth,  vers => $versToken,  
                      dir => "$auth/$versToken");
    my $structs = {
        G2R => {
            
        },
        R2P => {
        },
    };
    my $data = &gene_data();
    foreach my $gdat (values %{$data}) {
        my $gid = $gdat->{acc};
        my $gv  = $gdat->{v};
        foreach my $tdat (values %{$gdat->{rna}}) {
            foreach my $pdat (values %{$tdat->{prot}}) {
            }
        }
    }
}

sub go_ontology {
    
    my $data = &gene_data();

    my $nsj = "GeneOntology";
    my $type  = "Ontology";
    &go_meta();
    foreach my $nsi ("EnsemblGene", "EnsemblProtein") {
        my %fbits = (type => $type,  mod  => $specID,
                     ns1  => $nsi,   ns2  => $nsj,
                     auth => $auth,  vers => $versToken,  
                     dir => "$auth/$versToken");
        
        my $meta = {
            MatrixType => $fbits{type},
            Modifier   => $fbits{mod},
            Source     => $ftpDir,
            Namespace  => [$nsi, $nsj],
        };
        my $fmeta = { %{$stashMeta}, %{$meta} };
        my $trg   = &primary_path(%fbits);
        unless (&output_needs_creation($trg)) {
            &msg("Keeping existing $nsj file:", $trg);
            &post_process( %fbits, meta => $fmeta );
            next;
        }
        &msg("Structuring GeneOntology for $nsi");
        my (%goLinks, %objMeta);
        foreach my $gdat (values %{$data}) {
            my @dats = ($gdat);
            if ($nsi eq "EnsemblProtein") {
                @dats =();
                foreach my $tdat (values %{$gdat->{rna}}) {
                    push @dats, values %{$tdat->{prot}};
                }
            }
            foreach my $dat (@dats) {
                my $acc        = $dat->{acc};
                $goLinks{$acc} = [ map {$_, "UNK"} @{$dat->{link}{GO} || []} ];
                $objMeta{$acc} = {
                    Symbol      => $dat->{sym},
                    Description => $dat->{desc},
                };
            }
        }
        &make_go_matrix( \%goLinks, \%objMeta, $trg, \%fbits, $fmeta,
                         $nsUrl, $authLong, \&_citation_MTX, $taxDat,
            [qw(Symbol Description)], \%defColDef);
    }
}

sub gene_data {
    return $geneData if ($geneData);
    $globalMeta = {};
    my $tsv = &gene_file();
    &msg("Reading Gene metadata file", $tsv);
    $geneData = {};
    open(TSV, "<$tsv") || die "Failed to read Gene TSV file:\n  $tsv\n  $!";
    my $head = <TSV>;
    my ($lastGdat, $lastTdat, %counts);
    my $gb = $manifest->{build}[0] || "Unk";
    while (<TSV>) {
        s/[\n\r]+$//;
        my @row = split("\t");
        my ($gid, $sym, $coord, $tid, $pid, $desc, $vers) = @row;
        my $dat;
        if ($gid) {
            ## Gene row
            my $loc = "";
            if ($coord) {
                my ($chr, $s, $e, $str) = split(',', $coord);
                $loc = sprintf("%s.%s:%d..%d[%s]", $chr, $gb, $s, $e,
                               $str < 0 ? '-1' : '+1');
            }
            $dat = $lastGdat = $geneData->{$gid} = {
                acc   => $gid,
                sym   => $sym,
                desc  => $desc,
                coord => $loc,
                rna   => {},
                link  => {},
            };
            $globalMeta->{$gid} = {
                Symbol       => $dat->{sym},
                Description => $dat->{desc},
                Version     => $vers,
                Location    => $loc,
            };
            $counts{Gene}++;
        } elsif ($tid) {
            ## RNA row
            $dat = $lastTdat = $lastGdat->{rna}{$tid} = {
                acc   => $tid,
                sym   => $sym  || $lastGdat->{sym},
                desc  => $desc || $lastGdat->{desc},
                prot  => {},
                link  => {},
            };
            $globalMeta->{$tid} = {
                Symbol      => $dat->{sym},
                Description => $dat->{desc},
                Version     => $vers,
            };
            $counts{RNA}++;
        } elsif ($pid) {
            ## Protein row
            $dat = $lastTdat->{prot}{$pid} = {
                acc   => $pid,
                sym   => $sym  || $lastTdat->{sym},
                desc  => $desc || $lastTdat->{desc},
                link  => {},
            };
            $globalMeta->{$pid} = {
                Symbol      => $dat->{sym},
                Description => $dat->{desc},
                Version     => $vers,
            };
            $counts{Protein}++;
        } else {
            &err("Row lacking gene, RNA or protein ID");
        }
        ## Extract links
        for my $i (7..$#row) {
            if ($row[$i] =~ /^\/(\S+?)="(.+)"$/) {
                my ($k, $vals) = ($1, $2);
                $dat->{link}{$k} = [ split(',', $vals) ];
            }
        }
    }
    close TSV;
    &msg("Finished:", 
         "     Gene: $counts{Gene}",
         "      RNA: $counts{RNA}",
         "  Protein: $counts{Protein}");
    return $geneData;
}

sub gene_file {
    my $file = &local_file("GeneInfo-$ftpSpec-$vers.tsv");
    unless (&output_needs_creation($file)) {
        &msg("Using existing Gene file:", $file) unless ($tasks{MakeGene}++) ;
        return $file;
    }
    &msg("Creating Gene file for $specID $versToken");
    my $data = {
        genes => {},
        nowChr => "",
    };
    foreach my $fdat (@{$manifest->{files}}) {
        my ($emblFile, $chr) = @{$fdat};
        $data->{nowChr}  = $chr;
        &msg("  [<] $emblFile");
        my $fh   = IO::Uncompress::Gunzip->new( "$tmpDir/$emblFile" );
        my $feat = &_process_feat(); # Init with empty feature structure
        my $tag  = &_process_tag();  # Init with empty tag structure
        while (<$fh>) {
            # AFAICT all the useful gene data are in the feature table
            if (/^FT   (.+)[\n\r]+$/) {
                my $ft = $1;
                my $tagTxt = "";
                if ($ft =~ /^([a-z]\S+)\s+(.+)/i) {
                    ## New feature. Process old one, reset:
                    $feat = &_process_feat($feat, $data, $1, $2);
                } elsif ($ft =~ /^\s+\/([a-z_]+)="(.+)/i) {
                    ## New tag with quoted content
                    $tag = &_process_tag($tag, $feat, $1);
                    $tagTxt = $2;
                } elsif ($ft =~ /^\s+\/([a-z_]+)=(.+)/i) {
                    ## New tag with unquoted content. I sure hope
                    ## these are single lines...
                    $tag = &_process_tag({tag => $1, txt => [$2]}, $feat);
                } elsif ($ft =~ /^\s+\/.+=/) {
                    ## Uh. I seem to have overlooked a /key=val format
                    &msg("Unexpected key-val line: $ft");
                } elsif ($ft =~ /^\s+(\S.+)/) {
                    ## Extended tag text
                    $tagTxt = $1;
                }
                if ($tagTxt && $tag->{tag}) {
                    if ($tagTxt =~ /(.+)"$/) {
                        ## End of tag
                        push @{$tag->{txt}}, $1;
                        $tag = &_process_tag($tag, $feat);
                    } else {
                        ## Tag continues (hopefully) on next line
                        push @{$tag->{txt}}, $tagTxt;
                    }
                }
            }
        }
        &_process_feat($feat, $data);
        close $fh;
        ## warn "DEBUG!!"; last;
    }
    my @gids = sort keys %{$data->{genes}};
    &msg("  Done parsing EMBL files: ".scalar(@gids)." genes");
    open(OUT, ">$file") || die "Failed to write gene file:\n  $file\n  $!";
    my @head = qw(GeneID Symbol Coordinate Transcript Protein 
                  Description Version Links);
    print OUT join("\t", @head)."\n";
    foreach my $gid (@gids) {
        my $gdat = $data->{genes}{$gid};
        my @rids = sort keys %{$gdat->{rna} || {}};
        ## Cycle through the RNA data so we can assemble a list of protein IDs
        my (%pidH, @rpRows);
        foreach my $rid (@rids) {
            my $tdat = $gdat->{rna}{$rid};
            my @pids = sort keys %{$tdat->{prot} || {}};
            map { $pidH{$_} = 1 } @pids;
            my $rdesc =  &_unique_val($tdat->{desc}, "Desc for $rid") || "";
            my @rrow = ("","","", $rid, "", $rdesc,
                        $tdat->{v}, &_link_vals($tdat));
            #if ($rdesc =~ /transcript_id/) {
            #    die Dumper($tdat);
            #}
            
            ## RNA Metadata
            push @rpRows, \@rrow;
            &msg("[?] Multiple proteins for $rid: ".
                 join(' ', @pids)) if ($#pids > 0);
            foreach my $pid (@pids) {
                my $pdat = $tdat->{prot}{$pid};
                my @prow = ("","","", "", $pid,
                            &_unique_val($pdat->{desc}, "Desc for $pid"),
                            $pdat->{v}, &_link_vals($pdat));
                ## Protein Metadata
                push @rpRows, \@prow;
            }
        }
        my @pids = sort keys %pidH;
        my $gdesc =  &_unique_val($gdat->{desc}, "Desc for $gid");
        ## I've never found the parenthetical stuff at the end helpful, eg:
        ##   [Source:RGD Symbol;Acc:1308753]
        $gdesc =~ s/\s*\[[^\]]+\][\s\.]*$//;
        my @row = ($gid, 
                   &_unique_val($gdat->{sym}, "Symbol for $gid"), 
                   join(',', @{$gdat->{coord} || []}),
                   join(',', @rids),
                   join(',', @pids),
                   $gdesc,
                   $gdat->{v}, 
                   &_link_vals($gdat)
            );
        ## Gene Metadata
        print OUT join("\t", map { defined $_ ? $_ : "" } @row)."\n";
        ## print out the RNA and protein rows;
        foreach my $rpr (@rpRows) {
            print OUT join("\t", map { defined $_ ? $_ : "" } @{$rpr})."\n";
        }
    }
    close OUT;
    &msg("Generated Ensembl metadata TSV file", $file);
    return $file;
}

sub _link_vals {
    my $dat = shift;
    my @lk  = sort keys %{$dat->{link}};
    my @rv;
    foreach my $k (@lk) {
        my %u = map {$_ => 1} @{$dat->{link}{$k}};
        my @v = sort keys %u;
        push @rv, "/$k=\"".join(',', @v)."\"";
    }
    return @rv;
}

sub _unique_val {
    my ($arr, $context) = @_;
    my $num = 0;
    my %unq;
    map { $unq{$_} ||= ++$num } @{$arr};
    
    my @vals = sort {$unq{$a} <=> $unq{$b}} keys %unq;
    &msg("Non-unique values [$context]: ".join(' / ', @vals)) if ($num > 1);
    return $vals[0] || "";
}

sub _process_feat {
    my ($feat, $data, $newPtag, $coord) = @_;
    my $rv = {
        ptag  => $newPtag || "",
        coord => $coord, # We will only inspect for genes
        tags  => {},
    };
    if (my $fn = $feat->{ptag}) {
        ## This feature has some accumulated information that needs
        ## processing
        my (%desc, $sym);
        my $gCor = [ $data->{nowChr} ];
        my %tags = %{$feat->{tags} || {}};
        my @gids = @{$tags{gene} || []};
        my @tids = @{$tags{transcript_id} || []};
        my @pids = @{$tags{protein_id} || []};
        map { delete $tags{$_} } qw(gene transcript_id protein_id);
        if ($fn =~ /^(misc_RNA|mRNA)/) {
            push @tids, @{$tags{standard_name} || []};
            $desc{rna} = $tags{note} || [];
            map { delete $tags{$_} } qw(note standard_name);
        } elsif ($fn eq 'gene') {
            $sym = $tags{locus_tag} || [];
            $desc{gene} = $tags{note} || [];
            map { delete $tags{$_} } qw(note locus_tag);
            ## Let's also do some simple coordinate parsing
            if (my $coord = $feat->{coord}) {
                my $str = 1;
                if ($coord =~ /^complement\((.+)\)\s*$/) {
                    $coord = $1;
                    $str   = -1;
                }
                if ($coord =~ /^(\d+)\.\.(\d+)$/) {
                    $gCor = [ $gCor->[0], $1, $2, $str ];
                } else {
                    &msg("Unexpected gene coordinates: ".$feat->{coord});
                }
            }
        } elsif ($fn eq 'CDS') {
            foreach my $n (@{$tags{note}}) {
                if ($n =~ /^transcript_id=(.+)/) {
                    push @tids, $1;
                }
            }
            map { delete $tags{$_} } qw(note);
        } elsif (!$ignorePrimaryTag->{$fn} && !$tasks{UnkPriTag}{$fn}++) {
            &msg("Primary tag '$fn' does not have logic assigned to it");
        }

        ## Ok, we expect to at *least* have a gene ID if there's
        ## anyting useful to do. So start by seeing if we have one
        ## (and only one) ENSG accession.
        my %ug = map { $_ => 1 } @gids;
        @gids = keys %ug;
        ## Don't do anyting unless we have a gene ID:
        return $rv if ($#gids == -1);
        if ($#gids != 0) {
            &err("Multiple GeneIDs for $fn: ".join(' + ', @gids));
            return $rv;
        }
        my ($gid, $gv, $gt) = &_ok_ens_acc($gids[0]);
        if (!$gid) {
            &err("Unfamiliar gene ID format: ".$gids[0]);
            return $rv;
        }
        if ($gt ne 'G') {
            &err("GeneID does not have expected 'G': ".$gids[0]);
            return $rv;
        }
        my $gdat = $data->{genes}{$gid} ||= {
            acc   => $gid,
            coord => $gCor,
            rna   => {},
            sym   => [],
            link  => {},
        };
        push @{$gdat->{sym}}, @{$sym || []};
        $gdat->{desc} ||= $desc{gene};
        if ($gv) {
            ## Note the gene version, sanity check:
            if ($gdat->{v} && $gv != $gdat->{v}) {
                &err("GeneID with multiple versions: $gv + ".$gdat->{v});
            } else {
                $gdat->{v} ||= $gv;
            }
        }
        
        ## Do we have transcript and maybe protein IDs? If so, we can
        ## build the Gene-RNA-Protein hierarchy, and make some RNA- or
        ## protein-specific connections.
        my ($tdat, $pdat);
        if ($#tids > -1) {
            ## We have a transcript ID
            my %ut = map { $_ => 1 } @tids;
            @tids = keys %ut;
            if ($#tids > 0) {
                &err("Multiple Transcript IDs for $fn: ".join(' + ', @tids));
            } else {
                my ($tid, $tv, $tt) = &_ok_ens_acc($tids[0]);
                if (!$tid) {
                    &err("Unfamiliar transcript ID format: ".$tids[0]);
                } elsif ($tt ne 'T') {
                    &err("TranscriptID does not have expected 'T': ".$tids[0]);
                } else {
                    $tdat = $gdat->{rna}{$tid} ||= {
                        acc  => $tid,
                        prot => {},
                        link => {},
                    };
                    $tdat->{desc} ||= $desc{rna};
                    ## WEIRD. WHERE ARE THESE COMING FROM?? die "$tid = $desc{rna}" if ($desc{rna} && $desc{rna} =~ /^transcript_id/);
                    if ($tv) {
                        ## Note the RNA version, sanity check:
                        if ($tdat->{v} && $tv != $tdat->{v}) {
                            &err("Transcript with multiple versions: $tv + ".
                                 $tdat->{v});
                        } else {
                            $tdat->{v} ||= $tv;
                        }
                    }
                    ## Finally, can we associate a protein with the RNA?
                    my %up = map { $_ => 1 } @pids;
                    @pids = keys %up;
                    if ($#pids > 0) {
                         &err("Multiple Protein IDs for $fn: ".
                              join(' + ', @pids));
                    } elsif ($#pids == 0) {

                        
                        my ($pid, $pv, $pt) = &_ok_ens_acc($pids[0]);
                        if (!$pid) {
                            &err("Unfamiliar protein ID format: ".$pids[0]);
                        } elsif ($pt ne 'P') {
                            &err("Protein ID does not have expected 'P': ".
                                 $pids[0]);
                        } else {
                            $pdat = $tdat->{prot}{$pid} ||= {
                                acc  => $pid,
                                link => {},
                            };
                            if ($pv) {
                                ## Note the protein version, sanity check:
                                if ($pdat->{v} && $pv != $pdat->{v}) {
                                    &err("Protein with multiple versions: $pv + ".
                                         $pdat->{v});
                                } else {
                                    $pdat->{v} ||= $pv;
                                }
                            }
                        }
                    } # Protein substructure
                }
            } # Transcript substructure
        }
        
        foreach my $dx (@{$tags{db_xref}}) {
            if ($dx =~ /^([^:]+):(.+)/) {
                my ($db, $val) = ($1, $2);
                if ($db eq 'GO') {
                    if ($val =~ /^\d{7}$/) {
                        ## GO terms are generally associated with gene and prot:
                        push @{$gdat->{link}{GO}}, $dx;
                        push @{$pdat->{link}{GO}}, $dx if ($pdat);
                    } else {
                        &err("Weird GO accession: $dx");
                    }
                } elsif ($db =~ /^RefSeq_/) {
                    if ($db =~ /peptide/) {
                        push @{$pdat->{link}{RefSeqProtein}}, $val if ($pdat);
                    } else {
                        push @{$tdat->{link}{RefSeqRNA}}, $val if ($tdat);
                    }
                } elsif ($asIsXrefs->{$db}) {
                    if ($db =~ /^(MetaCyc|UniPathway)$/) {
                        ## Gene level association
                        push @{$gdat->{link}{$db}}, $val;
                    }
                    ## Associate all with protein, if a protein ID is known:
                    push @{$pdat->{link}{$db}}, $val if ($pdat);
                } elsif ($db eq 'Uniprot/SWISSPROT') {
                    push @{$pdat->{link}{UniProt}}, $val if ($pdat);
                    ## Note - we're ignoring TrEMBL
                } elsif ($db eq 'KEGG_Enzyme') {
                    push @{$gdat->{link}{KEGG}}, $val;
                    push @{$pdat->{link}{KEGG}}, $val if ($pdat);
                } elsif (!$ignoreXrefs->{$db}) {
                    my $msg = "Unknown DB $db in $fn";
                    &msg($msg) unless ($tasks{$msg}++);
                }
            }
        }
        delete $tags{db_xref};
        foreach my $tag (keys %tags) {
            unless ($ignoreFeatTags->{$tag}) {
                my $msg = "Un-utilized feature tag $tag in $fn";
                &msg($msg) unless ($tasks{$msg}++);
            }
        }
    }
    return $rv;
}

sub _ok_ens_acc {
    ## Check if an Ensembl accession is ok
    my $id = shift;
    if ($id =~ /^(ENS([A-Z]{3})?([GTP])\d{11})(\.(\d+))?$/) {
        ## Acc, Version, Type
        return ($1, $5, $3);
    }
    return "";
}

sub _process_tag {
    my ($tag, $feat, $newTag) = @_;
    if (my $tn = $tag->{tag}) {
        ## Structure the previously processed tag
        if (my $val = join(' ', @{$tag->{txt} || []})) {
            ## I've read the feature table specification:
            ##   http://www.insdc.org/files/feature_table.html
            
            ## I really can't see how one should know to join
            ## 'continuation lines' with or without a space. The one
            ## place I can see 'no space' being needed is following a
            ## dash, so I'll clean that up:
            $val =~ s/\- /\-/g;
            push @{$feat->{tags}{$tn}}, $val;
        }
    }
    return {
        tag => $newTag || "",
        txt => [],
    };
}

sub fetch_all_files {
    return undef if ($taxDat->{error});
    ## Ensembl breaks up EMBL files by chromosome.
    $ftpDir = sprintf("pub/release-%d/embl/%s/", $vers, $ftpSpec);
    ## Get the checksums file both as a check that we have a valid
    ## species and version, and as a manifest of tar.gz files:
    my $cksum = &fetch_url($ftpDir . "CHECKSUMS", "CHECKSUMS-$ftpSpec-$vers");
    unless (-s $cksum) {
        $taxDat->{error} ||= "Failed to find data at $ftpDir";
        return undef;
    }

    ## Parse checksum file
    open(CKS, "<$cksum") || die "Failed to read checksum file:\n  $cksum\n  $!";
    my $rv = {
        build   => {},
        species => {},
        files   => [],
    };
    while (<CKS>) {
        s/[\n\r]+$//;
        if (/^\d+\s\d+\s((.+?)\.(.+?)\.$vers\.(.+?).dat.gz)$/) {
            my ($file, $tax, $gb, $chr) = ($1, $2, $3, $4);
            $rv->{build}{$gb}++;
            $rv->{species}{$tax}++;
            if ($chr =~ /^chromosome\.(.+)$/) {
                $chr = $1;
            } elsif ($chr =~ /^(nonchromosomal)$/) {
                $chr = "OTHER";
            } else {
                warn "Unrecognized chromosome specifier '$chr'\n"
                    unless ($tasks{weirdChr}{$chr}++);
                $chr = "OTHER";
            }
            push @{$rv->{files}}, [$file, $chr];
            &fetch_url($ftpDir . $file);
        }
    }
    close CKS;
    ## Turn build and species count hashes into arrays. These *should*
    ## each yield a unique value.
    foreach my $key ('build', 'species') {
        my $h = $rv->{$key};
        $rv->{$key} = [ sort { $h->{$b} <=> $h->{$a} }  keys %{$h} ];
    }
    return $rv;
}

sub _citation_MTX {
    return &_default_parameter( "Citation", "Ensembl 2018; Zerbino et al; Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D754–D761; doi:10.1093/nar/gkx1098")."%\n";
}

