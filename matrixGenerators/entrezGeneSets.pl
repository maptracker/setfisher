#!/usr/bin/perl -w

use strict;
use LWP::UserAgent;
use Archive::Tar;
use File::Path 'mkpath';
use IO::Uncompress::Gunzip;

my $defTmp = "/tmp/entrezGeneSets";
my $defFtp = "https://ftp.ncbi.nih.gov/";
my $defEC  = "P,IEA,NAS,IRD,IBD,IBA,RCA,IGC,ISS,ISA,ISM,ISO,TAS,EXP,IEP,IPI,IMP,IGI,IDA";

my $args = &parseargs({
    tmpdir   => $defTmp,
    ftp      => $defFtp,
    eclevels => $defEC,
    dir      => ".",
    clobber  => 0,
    verify_hostname => 0,
});

my $geneIdUrl  = 'https://www.ncbi.nlm.nih.gov/gene/%s'; # For integer IDs
my $symUrl     = 'https://www.ncbi.nlm.nih.gov/gene/?term=%s%%5Bsym%%5D';
my $locLinkUrl = 'https://www.ncbi.nlm.nih.gov/gene/?term=%s'; # For LOC###
my $bar        = "- " x 20;
my $mtxSep     = " :: ";

=head2 SSL Issues with NCBI

At time of writing (9 May 2017) there was an issue with a mismatch
between the FTP site and the certificate:

   ERROR: certificate common name “*.ncbi.nlm.nih.gov” doesn’t match requested host name “ftp.ncbi.nih.gov”.

This will cause LWP to fail unless verify_hostname is set to false,
which is the current default. When NCBI supplies a certificate with
matching common name this option should be removed.

=cut

my $tmpDir   = $args->{tmpdir}; $tmpDir =~ s/\/+$//;
my $outDir   = $args->{dir};    $outDir =~ s/\/+$//;
my $ftpUrl   = $args->{ftp};    $ftpUrl =~ s/\/+$//;
my $clobber  = $args->{clobber} || 0;
my $ua       = LWP::UserAgent->new;
$ua->env_proxy;
if (my $prox = $args->{proxy}) {
    $ua->proxy(['http', 'ftp'], $prox);
    &msg("Set proxy:", $prox);
}
foreach my $sslopt (qw(verify_hostname SSL_ca_file SSL_ca_path)) {
    ## Certificates are fun!
    my $v = $args->{lc($sslopt)};
    if (defined $v) {
        $ua->ssl_opts( $sslopt, $v);
        &msg("SSL Option:", "$sslopt => $v");
    }
}

&mkpath([$tmpDir, $outDir]);

my $taxDat   = &extract_taxa_info( $args->{species} );

if ($taxDat->{error} || $args->{h} || $args->{help}) {
    warn "
Usage:

This program will generate gene sets for the SetFisher enrichment
package. It will recover information from the Entrez FTP site based on
a species identifier you provide.

Optional Arguments:

      -ftp Default '$defFtp'. URL to Entrez's FTP site

   -tmpdir Default '$defTmp'. Directory holding downloaded files

      -dir Default '.'.  Output directory that will contain generated
           matrix files

  -clobber Default 0, which will preserve any already-generated
           files. Set to 1 to regenerate matrix files, and any value
           greater than 1 to also re-download base files from Entrez.

     -name A name to use for the files. By default this will be chosen
           from available taxonomy information, preferring the Genbank
           Common Name

 -eclevels Default '$defEC'. This series maps evidence codes to an
           integer 'score' that will be used to represent each
           locus->GO assignment in the MatrixMarket file. 'Worse'
           evidence codes should be listed first. The default ranking
           is generally reasonable, but also fairly arbitrary on a
           fine scale.

";
    warn "\n[!!] Failed to find target species:\n     $taxDat->{error}\n\n"
        if ($taxDat->{error});
    exit;
}
my ($geneMeta);

my $taxid = $taxDat->{taxid};


&msg("Working directory:", $tmpDir);
&msg("Output directory:", $outDir);

## Summarize the species we are going to parse
my @taxBits = ("TaxID: $taxid");
my $specID  = $args->{name}; # This will be a file naming stub
if (my $gcn = $taxDat->{'genbank common name'}) { 
    push @taxBits, "Common Name: ".join(' // ', @{$gcn});
    unless ($specID) {
        ## In most cases, the GCN will be the file token:
        $specID = $gcn->[0];
        ## dog -> Dog
        substr($specID, 0, 1, uc( substr($specID, 0, 1) ));
    }
    
}
if (my $sn = $taxDat->{'scientific name'}) { 
    push @taxBits, "Scientific Name: ".join(' // ', @{$sn});
    $specID ||= $sn->[0];
}
if (my $cn = $taxDat->{'common name'}) { 
    push @taxBits, "Other Aliases: ".join(' // ', @{$cn});
}
$specID ||= "taxa$taxid";
$specID =~ s/\s+/_/g;
push @taxBits, "File token: $specID";

&msg("Target species:", @taxBits);

&map_symbol();
&ontology_go();
&ontology_pubmed();

sub make_metadata_file {
    # Simple TSV file of species-specific metadata
    my $trg = sprintf("%s/Metadata-%s_Entrez.tsv", $outDir, $specID);

    unless (&output_needs_creation($trg)) {
        &msg("Using existing Metadata file:", $trg);
        return $trg;
    }
    open (METAF, ">$trg") || &death("Failed to wirte metadata file",
                                    $trg, $!);
    
    &msg("Parsing basic Gene metadata");
    ## If you see 'expected column' errors, you will need to change
    ## the right hand value (after the '=>') of the offending column,
    ## after determining what the new column name is:
    my $cols = { 
        TaxID       => 'tax_id',
        GeneID      => 'GeneID',
        Symbol      => 'Symbol',
        AuthSym     => 'Symbol_from_nomenclature_authority',
        Type        => 'type_of_gene',
        Status      => 'Nomenclature_status',
        Aliases     => 'Synonyms',
        Description => 'description',
    };
    my @colOut = qw(GeneID Symbol Type Description Aliases);
    print METAF join("\t", @colOut, "SymScore") ."\n";

    my ($fh)  = &gzfh("gene/DATA/gene_info.gz", "gene_info.gz", $cols);
    my @getInds = map { $cols->{$_} } @colOut;
    my $tInd = $cols->{TaxID};   # Taxid, eg 9606
    my $lInd = $cols->{GeneID};  # GeneID, eg 859
    my $sInd = $cols->{Symbol};  # The 'main' symbol, eg CAV3
    my $aInd = $cols->{Aliases}; # Alt symbols, eg LGMD1C|LQT9|VIP-21|VIP21
    my $oInd = $cols->{Status};  # Nom.status, eg O, I or blank
    my $nInd = $cols->{AuthSym}; # Nom.auth symbol
    my $ngene = 0;

    while (<$fh>) {
        s/[\n\r]+$//;
        my @row = split("\t");
        unless ($row[$tInd] == $taxid) {
            ## Stop parsing if we have already found some hits. THIS
            ## MAY BE A BAD IDEA. It should speed up parsing (neglects
            ## need to scan whole file) but there is no guarantee that
            ## taxa won't end up scattered through file in future
            ## versions:
            last if ($ngene);
            next; # Otherwise, keep scanning for taxid
        }
        map { s/^\-$// } @row; # Turn '-' cells to empty string
        my $gid = $row[ $lInd ];
        unless ($gid) {
            &err("Row without GeneID ???", $_);
            next;
        }
        # I am somewhat unclear on the symbol situation. The 'Symbol'
        # column is NOT unique (eg RNR1 and RNR2 in human). The
        # Symbol_from_nomenclature_authority is not always reflected
        # in Symbol or Synonyms (eg MT-RNR1 is listed as RNR1). So I
        # am going to collect all symbols uniquely in Aliases, and set
        # the Symbol to be:
        #    Symbol_from_nomenclature_authority
        #    Symbol
        #    First Alias

        my $aliText  = $row[ $aInd ] || "";
        # Symbols (at least aliases) *can* contain whitespace. Example:
        #   https://www.ncbi.nlm.nih.gov/gene/137964 -> "1-AGPAT 6"
        # ~168 human aliases have a space in them.
        # So just remove whitespace on edges of aliases:
        $aliText     =~ s/^\s+//; $aliText =~ s/\s+$//; # Edge whitespace
        # Found at least one instance where a stray semicolon snuck in
        #   https://www.ncbi.nlm.nih.gov/gene/7334 -> "UBCHBEN; UBC13"
        # ... and one with a comma:
        #   https://www.ncbi.nlm.nih.gov/gene/441282 -> ""
        # Pretty sure the semicolon is not allowed in a symbol, but it
        # looks like the comma can be:
        #   https://www.ncbi.nlm.nih.gov/gene/10317 -> "beta-1,3-GalTase 5"
        $aliText     =~ s/\s*[;]\s*/\|/g;
        # I don't think the following is a problem (eg 'ABC||XYZ') but
        # will grep it out just in case:
        $aliText     =~ s/\|[\s\|]+/\|/g; # Null entries
        $aliText     = "" if ($aliText eq '|'); # Nothing at all
        # Break into individual symbols:
        my @aliases  = split(/\s*\|\s*/, $aliText);
        my %uAli     = map { $aliases[$_] => $_ + 1 } (0..$#aliases);

        # The 'main' "Symbol" column:
        my $sym      = $row[ $sInd ] || "";
        my $stat     = $row[ $oInd ] || "";
        my $priScore = $stat eq "O" ? "1.0" : $stat eq "I" ? "0.5" : "0.4";

        if (my $nas = $row[ $nInd ]) {
            ## Nomenclature Authority (eg HGNC) symbol is present
            if ($nas ne $sym) {
                ## The authoritative symbol is different that the main one
                ## https://www.ncbi.nlm.nih.gov/gene/4549
                ## Official Sym = MT-RNR1, gene_info.gz Sym = RNR1
                ## ... which collides with:
                ## https://www.ncbi.nlm.nih.gov/gene/6052
                ## Official Sym = gene_info.gz Sym = RNR1
                ## Make sure the Symbol column is in the aliases:
                $uAli{$sym} ||= 0;
                ## And reset the displayed Symbol:
                $sym = $row[ $sInd ] = $nas;
            }
        }
        unless ($sym) {
            ## No symbol set
            if (my $usym = $aliases[0]) {
                ## There are aliases - take the first one as "the" symbol:
                $sym      = $usym;
                $priScore = 0.4;
                delete $uAli{$sym};
            }
        }
        delete $uAli{""};
        delete $uAli{$sym}; # Do not include the "main" symbol in the aliases
        # Reassemble the aliases in original order:
        $row[ $aInd ] = join
            ('|', sort { $uAli{$a} <=> $uAli{$b} } keys %uAli) || "";
        
        ## Add the entry to the TSV file:
        print METAF join("\t", (map { $row[$_] } @getInds), $priScore)."\n";
        $ngene++;
    }
    close $fh;
    close METAF;
    &msg("Generated Entrez metadata file", $trg);
    return $trg;
}

sub capture_gene_meta {
    return $geneMeta if ($geneMeta);
    $geneMeta  = {};
    my $mdFile = &make_metadata_file();
    open(METAF, "<$mdFile") || &death("Failed to read metadata TSV",
                                      $mdFile, $!);
    my $head = <METAF>;
    $head =~ s/[\n\r]+$//;
    my @cols = split(/\t/, $head);
    while (<METAF>) {
        s/[\n\r]+$//;
        my @row = split(/\t/);
        my %data = map { $cols[$_] => $row[$_] } (0..$#cols);
        my $gid  = $data{GeneID};

        if ($geneMeta->{$gid}) {
            &err("Multiple entries for GeneID == $gid");
        } else {
            ## Turn alias string into array
            my @aliases = split(/\|/, $data{Aliases} || "");
            @aliases = () if ($#aliases == 0 && $aliases[0] eq '');
            $data{Aliases} = \@aliases;
            $geneMeta->{$gid} = \%data;
        }
    }
    close METAF;
    return $geneMeta;
}

sub map_symbol {
    my $trg = sprintf("%s/Map-%s_Symbol-to-Entrez.mtx", 
                      $outDir, $specID);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing GeneOntology file:", $trg);
        return $trg;
    }
    &capture_gene_meta();
    my (%symbols, %genes);
    my ($gnum, $nznum) = (0, 0);
    my $scoreTags = {
        Official => "1.0",
        Interim  => "0.5",
        Unofficial => "0.4 (preferred) or 0.3",
    };

    # Symbols should all be the same case - but I'm not certain of
    # that. They will be collected under an upper-case key, but
    # preserve the first observed case for recording in the matrix
    my (%statuses, %oddChar);
    foreach my $gid (sort { $a <=> $b } keys %{$geneMeta}) {
        my $meta = $geneMeta->{$gid};
        my $pri  = $meta->{Symbol} || "";
        my $ns = 0;
        my $gn;
        if ($pri) {
            # "Main" / "primary" symbol
            my $pSc  = ($meta->{SymScore}  || 0.4) + 0;
            my $targ = $symbols{uc($pri)} ||= {
                name => $pri,
                hits => [],
            };
            $gn ||= $genes{$gid} ||= ++$gnum;
            $statuses{ $pSc == 1 ? "Official" :
                         $pSc == 0.5 ? "Interim" : "Unofficial" }++;
            push @{$targ->{hits}}, ($gn, $pSc);
            $nznum++;
            $ns++;
        }
        my @alis = @{$meta->{Aliases}};
        my $nUn  = $#alis + 1;
        $ns     += $nUn;
        $statuses{"Unofficial"} += $nUn;

        for my $i (0..$#alis) {
            ## All aliases will have a score of 0.3
            my $sym = $alis[$i];
            my $targ = $symbols{uc($sym)} ||= {
                name => $sym,
                hits => [],
            };
            $gn ||= $genes{$gid} ||= ++$gnum;
            push @{$targ->{hits}}, ($gn, 0.3);
            $nznum++;
            ## Tally odd characters. Mostly curiosity, but these may
            ## cause issues in some workflows.
            ## Allowed atypical characters:
            ##   '_' eg C4B_2 -> https://www.ncbi.nlm.nih.gov/gene/100293534
            ##   '@' eg HOXA@ -> https://www.ncbi.nlm.nih.gov/gene/3197
            $sym =~ s/[a-z0-9\.\-_@]//gi;
            foreach my $char (split('', $sym)) {
                if ($char) {
                    $oddChar{$char}{n}++;
                    $oddChar{$char}{ex} ||= $alis[$i];
                }
            }
        }
    }
    my @counts;
    foreach my $sdat (values %symbols) {
        my $ngene = ($#{$sdat->{hits}} + 1) / 2;
        $counts[$ngene < 10 ? $ngene : 10 ]++;
    }

    my @rowIds  = sort keys %symbols;
    my @gids    = sort { $genes{$a} <=> $genes{$b} } keys %genes;
    my ($rnum,$cnum) = ($#rowIds + 1, $#gids + 1);

    open(MTX, ">$trg") || &death("Failed to write Symbol map", $trg, $!);
    print MTX "%%MatrixMarket matrix coordinate real general
% Mapping from $specID symbols to Entrez IDs
% $rnum x $cnum sparse matrix with $nznum non-zero cells
% -- The 'setfisher' package can parse these comments to decorate the matrix --
% Separator '$mtxSep'
%% DEFAULT Name $specID Symbol-to-Entrez Map
%% DEFAULT Description 1.0=Official Symbol, 0.5=Interim, 0.4=Preferred unofficial, 0.3=Unofficial
% 'Preferred unofficial' = Just the first listed unofficial symbol
%
% Rows are Gene Symbols
%% DEFAULT RowDim Gene Symbol
%% DEFAULT RowUrl $symUrl
% Columns are Entrez IDs
%% DEFAULT ColDim Entrez ID
%% DEFAULT ColUrl $geneIdUrl
%
% CAUTION: Some symbols have case/capitalization subtleties. While
% there are often historical reasons for these case decisions,
% attempting to extract information from a gene symbol's case is about
% as reliable as determining the contents of a wrapped present by
% listening to the noises it makes when shaken. If you are pivotting
% data from symbols (rather than accessions) you're already working at
% a significant disadvantage:
%             PLEASE MATCH SYMBOLS CASE-INSENSITIVELY
% Case is preserved to aid in 'pretty' display of the symbols.
%
% $bar
% Symbol statuses:
";

    # Basic statistics commentary
    foreach my $stat (sort { $statuses{$b} <=> $statuses{$a} } keys %statuses) {
        printf(MTX "%% %15s : %d : score %s\n", $stat, $statuses{$stat},
               $scoreTags->{$stat} || "??");
    }
    print MTX "%
% Number of genes refererenced by symbol := Number of symbols
";
    for my $i (1..$#counts) {
        printf(MTX "%%  %s := %d\n", $i == 10 ? ">9" : " $i", $counts[$i] || 0);
    }
    my @ocs = sort { $oddChar{$b}{n} <=> $oddChar{$a}{n} } keys %oddChar;
    if ($#ocs != -1) {
        print MTX "%
% Potentially troublesome character := Number of symbols (Example)
";
        foreach my $oc (@ocs) {
            my $ocd = $oddChar{$oc};
            printf(MTX "%%  '%s' := %5d (%s)\n", $oc, $ocd->{n}, $ocd->{ex});
        }
    }
    print MTX "% $bar
% Comment blocks for [Row Name] and [Col Name] follow, followed finally by
% the triples that store the actual mappings.
% $bar
";

    ## Note row metadata
    my @meta = qw(Official NotOfficial);
    printf(MTX "%% Row %s\n", join($mtxSep, "Name", @meta));
    for my $i (0..$#rowIds) {
        my $dat  = $symbols{$rowIds[$i]};
        my $hits = $dat->{hits};
        my ($o, $u) = (0,0);
        for (my $j = 1; $j <= $#{$hits}; $j += 2) {
            if ($hits->[$j] < 1) { $u++ } else { $o++ }
        }
        printf(MTX "%% %d %s\n", $i+1, join($mtxSep, $dat->{name}, $o, $u));
    }

    ## Column metadata
    print MTX "\% $bar\n";
    &mtx_entrez(\@gids, 'Col');

    ## Finally, note triples for non-zero cells
    print MTX "% $bar
% Matrix triples : Row Col Score
% $bar
";
    printf(MTX "  %d %d %d\n", $rnum, $cnum, $nznum);
    for my $i (0..$#rowIds) {
        my $hits = $symbols{$rowIds[$i]}{hits};
        for (my $j = 0; $j < $#{$hits}; $j += 2) {
            printf(MTX "%d %d %0.1f\n", $i+1, $hits->[$j], $hits->[$j+1]);
        }
    }
    close MTX;
    &msg("Generated Symbol to Entrez mapping", $trg);
    return $trg;
}

sub mtx_entrez {
    my ($ids, $rc) = @_;
    &capture_gene_meta();
    my @meta = qw(Symbol Type Description);
    printf(MTX "%% %s %s\n", $rc, join($mtxSep, "Name", @meta));
    for my $i (0..$#{$ids}) {
        my $id = $ids->[$i];
        my $m = $geneMeta->{$id} || {};

        my @line = ($id, map { defined $_ ? $_ : "" } map { $m->{$_} } @meta);
        printf(MTX "%% %d %s\n", $i+1, join($mtxSep, @line));
    }
}

sub ontology_pubmed {
    my $trg = sprintf("%s/Ontology-%s_Entrez-to-PubMed.mtx", 
                      $outDir, $specID);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing GeneOntology file:", $trg);
        return $trg;
    }
    &capture_gene_meta();
    &msg("Parsing PubMed");
    ## If you see 'expected column' errors, you will need to change
    ## the right hand value (after the '=>') of the offending column,
    ## after determining what the new column name is:
    my $cols = { 
        taxid  => 'tax_id',
        gene   => 'GeneID',
    };
    my ($fh) = &gzfh("gene/DATA/gene2pubmed.gz", "gene2pubmed.gz", $cols);
    close $fh;
    warn "     ... writing matrix market file ...\n";
    return $trg;
}

sub ontology_go {
    my $trg = sprintf("%s/Ontology-%s_Entrez-to-GeneOntology.mtx", 
                      $outDir, $specID);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing GeneOntology file:", $trg);
        return $trg;
    }
    &capture_gene_meta();
    &msg("Parsing GeneOntology");
    ## If you see 'expected column' errors, you will need to change
    ## the right hand value (after the '=>') of the offending column,
    ## after determining what the new column name is:
    my $cols = { 
        taxid  => 'tax_id',
        gene   => 'GeneID',
        goid   => 'GO_ID',
        ec     => 'Evidence',
        term   => 'GO_term',
        ## pmid   => 'PubMed', # Can not really record this in MTX
        cat    => 'Category',
    };

    ## Set up EvidenceCode -> 'score' map
    my $ecTxt = uc($args->{eclevels} || "");
    $ecTxt =~ s/^\s+//; $ecTxt =~ s/\s+$//; # leading/trailing whitespace
    my @ecl = split(/\s*,\s*/, $args->{eclevels});
    my %ecSc = map { $ecl[$_] => $_ + 1 } (0..$#ecl);

    my ($fh) = &gzfh("gene/DATA/gene2go.gz", "gene2go.gz", $cols);
    my $tInd = $cols->{taxid}; # Taxid, eg 9606
    my $gInd = $cols->{goid};  # GO ID, eg GO:0005886
    my $eInd = $cols->{ec};    # Evidence code, eg TAS
    my $lInd = $cols->{gene};  # Entrez Gene ID, eg 859
    my (%genes, %goMeta);
    my ($nznum, $ngo) = (0,0);

    while (<$fh>) {
        s/[\n\r]+$//;
        my @row = split("\t");
        unless ($row[$tInd] == $taxid) {
            ## Stop parsing if we have already found some hits. THIS
            ## MAY BE A BAD IDEA. It should speed up parsing (neglects
            ## need to scan whole file) but there is no guarantee that
            ## taxa won't end up scattered through file in future
            ## versions:
            last if ($ngo);
            next; # Otherwise, keep scanning for taxid
        }
        map { s/^\-$// } @row; # Turn '-' cells to empty string
        my $gid = $row[$gInd];
        ## Capture term metadata once
        my $gm = $goMeta{$gid} ||= {
            Description => $row[ $cols->{term} ],
            Category    => $row[ $cols->{cat}  ],
            order       => ++$ngo,
        };
        my $sc = $ecSc{ $row[ $eInd ] };
        unless ($sc) {
            unless (defined $sc) {
                &err("Unknown Evidence code '$row[ $eInd ]'",
                     "This code was not in -eclevels, it will be ignored");
                $ecSc{ $row[ $eInd ] } = 0;
            }
            next;
        }
        # Capture orderID/score pairs
        push @{$genes{ $row[ $lInd ] }}, ( $gm->{order}, $sc );
        $nznum++;
    }
    close $fh;
    warn "     ... writing matrix market file ...\n";

    my @goids = sort { $goMeta{$a}{order} <=> $goMeta{$b}{order} } keys %goMeta;
    my @gids  = sort { $a <=> $b } keys %genes;
    my ($rnum, $cnum) = ($#gids + 1, $#goids + 1);

    open(MTX, ">$trg") || &death("Failed to write Symbol map", $trg, $!);
    print MTX "%%MatrixMarket matrix coordinate real general
% Mapping from $specID Entrez IDs to GeneOntology Terms
% $rnum x $cnum sparse matrix with $nznum non-zero cells
% -- The 'setfisher' package can parse these comments to decorate the matrix --
% Separator '$mtxSep'
%% DEFAULT Name $specID Entrez GeneOntology
%% DEFAULT Description GeneOntology assignments for $specID Entrez genes
%
% Rows are Entrez IDs
%% DEFAULT RowDim Entrez ID
%% DEFAULT RowUrl $geneIdUrl
% Columns are GO Terms
%% DEFAULT ColDim GO Term
%% DEFAULT ColUrl http://amigo.geneontology.org/amigo/term/%s
%
% Terms with few genes assigned to them struggle to reach statistical
% significance. Excluding them removes some distractions, but also
% helps minimize multiple-testing penalties on 'long shot' terms. The
% threshold below represents the minimum number of genes to be
% assigned to a term for the term to be kept.
%
%% DEFAULT minSetSize 7
%
% A major distortion in Fisher-based testing can come from
% 'questionable' genes. These include real-but-untranscribed entities
% (eg pseudogenes), speculative entries that will eventually be
% removed from the transcriptome, or rare genes that are expressed
% only in certain tissues or developmental stages. From a
% marbles-in-an-urn perspective, all these categories represent
% marbles that can NEVER be removed from the urn - the effect is to
% inflate the significance of those that you do, since the world size
% is larger than it really should be. A fairly simple way to
% automatically recognize such objects is to look at the total number
% of terms annotated to each gene, since speculative genes tend to be
% unannotated. The filter below removes genes with fewer than the
% indicated number of terms assigned to it.
%
%% DEFAULT minOntoSize 2
%
% High-level GO terms tend to be less useful in biological
% interpretation. We have also observed that they tend to be
% disproportionately significant as they are greatly impacted by the
% effect mentioned above (un-selectable genes distorting the
% significance of enriched sets). Terms that are assigned to a higher
% percentage of the world below will be excluded from the analysis. In
% addition to possibly bringing some clarity to reports (how often is
% a significant hit to 'Catalytic process' going to help you?) it
% brings a small multiple testing benefit.
%
%% DEFAULT maxSetPerc 5
%
% Used in SetFisher, the filters above would be applied recursively;
% Terms are removed, then genes, and the process is repeated until no
% further alterations occur.
%
%% $bar
%% Values should be treated as factors representing evidence codes:
";

    my $top = "%% ------     ";
    my $bot = "%% LEVELS [,][";
    for my $li (0..$#ecl) {
        my $l = $ecl[$li];
        my $v = $ecSc{$l};
        $l .= ',' unless ($li == $#ecl);
        $bot .= $l;
        $top .= sprintf('%-'.CORE::length($l).'s', $v);
    }
    print MTX "$top\n";
    print MTX "$bot]\n";
    print MTX "%% Lower factor values should represent lower confidence evidence.\n";
    print MTX "%% http://geneontology.org/page/guide-go-evidence-codes\n";
    print MTX "% $bar
% Comment blocks for [Row Name] and [Col Name] follow, followed finally by
% the triples that store the actual mappings.
% $bar
";

    ## Row metadata
    &mtx_entrez(\@gids, 'Row');
    print MTX "\% $bar\n";

    ## Column metadata
    my @meta = qw(Description Category);
    printf(MTX "%% Col %s\n", join($mtxSep, "Name", @meta));
    for my $i (0..$#goids) {
        my $goid = $goids[$i];
        my $dat  = $goMeta{$goid};
        printf(MTX "%% %d %s\n", $i+1, join($mtxSep, $goid,
                                            map { $dat->{$_} } @meta));
    }

    ## Finally, note triples for non-zero cells
    print MTX "% $bar
% Matrix triples : Row Col Score
% $bar
";
    printf(MTX "  %d %d %d\n", $rnum, $cnum, $nznum);
    for my $i (0..$#gids) {
        my $hits = $genes{$gids[$i]};
        for (my $j = 0; $j < $#{$hits}; $j += 2) {
            printf(MTX "%d %d %d\n", $i+1, $hits->[$j], $hits->[$j+1]);
        }
    }
    close MTX;
    &msg("Generated Entrez GO ontology", $trg);

    die "Working here";
    return $trg;
}

sub parseargs {
    ## Command line argument parsing
    my $rv = $_[0] || {};
    my $i = 0;
    while ($i <= $#ARGV) {
        my $key = lc($ARGV[$i]);
        $key =~ s/^\-+//;
        my $val = $i < $#ARGV && $ARGV[$i+1] !~ /^\-/ ? $ARGV[++$i] : 1;
        $rv->{$key} = $val;
        $i++;
    }
    return $rv;
}

sub msg {
    warn "[*] ".join("\n    ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}

sub err {
    warn "[!!] ERROR: ".join
        ("\n     ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}

sub death { &err(@_); die " -- "; }

sub gzfh {
    my ($urlDir, $locFile, $expectCols) = @_;
    my $src = &fetch_url($urlDir, $locFile);
    my $fh  = IO::Uncompress::Gunzip->new( $src ) ||
        &death("Failed to gunzip file", $src);
    if ($expectCols) {
        # We are expecting particular columns to be present
        # Verify and set the column index when we find them
        my $head = <$fh>;
        $head =~ s/^#+//; # Initial number sign on most headers
        $head =~ s/[\n\r]+$//;
        my @found = split(/\t/, $head);
        my %lu = map { $found[$_] => $_ } (0..$#found);
        while (my ($tok, $col) = each %{$expectCols}) {
            my $ind = $lu{$col};
            if (defined $ind) {
                $expectCols->{$tok} = $ind;
            } else {
               &death("Failed to find expected column: '$col'",
                      "The file format may have changed. Please check it:",
                      $src, "... and then scan this program for '\$cols'",
                      "That variable defines expected columns after each '=>'",
                      "Find the offending value, and replace it with the value observed in the .gz file",
                      "Bear in mind there are several subroutines with '\$cols' - find the relevant one.");
            }
        }
    }
    return $fh;
}

sub extract_taxa_info {
    ## Used to deconvolute user species request into a formal
    ## species. In particular, we'll need the taxid (eg 9606 for
    ## human) to extract relevant subsets of information.
    my $req = shift;
    return { error => "No taxa request specified" }  unless ($req);
    my $srcFile = "$tmpDir/names.dmp";
    if (&source_needs_recovery($srcFile)) {
        my $tgz = &fetch_url("pub/taxonomy/taxdump.tar.gz", "taxdump.tar.gz");
        my $tar = Archive::Tar->new($tgz);
        $tar->extract_file("names.dmp", $srcFile);
        die join("\n  ", "Failed to extract taxa names",
                 $srcFile, "") unless (-s $srcFile);
        &msg("Extracted Taxonomy names", $srcFile);
    }
    open(TAX, "<$srcFile") || die 
        join("\n", "Failed to parse taxonomy file", $srcFile, $!, "");
    my $obj = { taxid => 0 };

    while (<TAX>) {
        my ($txid, $name, $uniq, $cls) = split(/\s*\|\s*/);
        if ($txid != $obj->{taxid}) {
            ## New name block
            if ($obj->{found}) {
                ## This is what we came for
                close TAX;
                return $obj;
            }
            $obj = { taxid => $txid };
            $obj->{found} = "tax_id" if ($req eq $txid);
        }
        push @{$obj->{$cls}}, $name;
        $obj->{found} = $cls if (lc($req) eq lc($name));
    }
    close TAX;
    return { error => "No match for '$req' found in $srcFile" };
}

sub fetch_url {
    my ($uReq, $dReq) = @_;
    my $dest = "$tmpDir/$dReq"; # Local file path
    if (&source_needs_recovery($dest)) {
        ## File not yet downloaded, or request to re-download
        my $url = "$ftpUrl/$uReq"; # Remote URL
        my $res = $ua->get( $url, ':content_file' => $dest );
        if (-s $dest) {
            &msg("Downloaded $dReq", $dest);
        } else {
            die join("\n  ", "Faliled to recover file",
                     "Source: $url", "Destination: $dest",
                     sprintf("HTTP Result: %s=%s",
                             $res->code(), $res->message()),  "");
        }
    }
    return $dest;
}


sub output_needs_creation {
    my $path = shift;
    ## Not if the file exists, is non-zero size, and clobber is false
    return 0 if (!$clobber && (-s $path));
    ## Otherwise yes:
    return 1;
}

sub source_needs_recovery {
    my $path = shift;
    ## Not if the file exists, is non-zero size, and clobber is 0 or 1
    return 0 if ($clobber < 2 && (-s $path));
    ## Otherwise yes
    return 1;
}
