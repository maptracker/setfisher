#!/usr/bin/perl -w

use strict;
my $scriptDir;
our $defTmp = "/tmp/entrezGeneSets";
my $defFtp = "ftp.ncbi.nih.gov";
my $defEC  = "P,IEA,NAS,IRD,IBD,IBA,RCA,IGC,ISS,ISA,ISM,ISO,TAS,EXP,IEP,IPI,IMP,IGI,IDA";
my $defSym = "Unknown,Unofficial,UnofficialPreferred,Interim,Official";

## RefSeq Status codes:
## https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_status_codes/
my $defRsStat = "SUPPRESSED,NA,WGS,MODEL,INFERRED,PREDICTED,PROVISIONAL,REVIEWED,VALIDATED";


our $defaultArgs = {
    ftp       => $defFtp,
    eclevels  => $defEC,
    rslevels  => $defRsStat,
    symlevels => $defSym,
    dir       => ".",
};

BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}
use lib "$scriptDir";
require Utils;
our ($args, $clobber, $ftp, $tmpDir, $maxAbst);
my ($dbh, $getPMID, $bfdbh, $bfClear, $bfSet);

use DBI;
use DBD::SQLite;
use Archive::Tar;
use XML::Twig;

my $geneIdUrl  = 'https://www.ncbi.nlm.nih.gov/gene/%s'; # For integer IDs
my $symUrl     = 'https://www.ncbi.nlm.nih.gov/gene/?term=%s%%5Bsym%%5D';
my $locLinkUrl = 'https://www.ncbi.nlm.nih.gov/gene/?term=%s'; # For LOC###
my $bar        = "- " x 20;

my $outDir   = $args->{dir};    $outDir =~ s/\/+$//;

&mkpath([$outDir]);

my $taxDat   = &extract_taxa_info( $args->{species} || $args->{taxa} );

if ($taxDat->{error} || $args->{h} || $args->{help}) {
    warn "
Usage:

This program will generate gene sets for the SetFisher enrichment
package. It will recover information from the Entrez FTP site based on
a species identifier you provide.

Required Argument:

  -species The species name you wish to extract

Optional Arguments:

      -ftp Default '$defFtp'. URL to Entrez's FTP site

   -tmpdir Default '$defTmp'. Directory holding downloaded files

      -dir Default '.'.  Output directory that will contain generated
           matrix and metadata files

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

 -rslevels Default '$defRsStat'.
           Like -eclevels, but maps RefSeq status codes to an integer
           score. Small values are considered 'lower confidence'.

-symlevels Default '$defSym'.
           Factor levels for symbol nomenclature status.

           See: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_status_codes/

 -pubmeddb A path to a SQLite database that holds PubMed (Medline)
           publication dates and article titles. While optional, the
           PubMed ontology will not be created in the absence of this
           database (it's really not too informative to only know the
           PubMed ID for an enriched publication).

           To create this supporting database, see the included utility:

                 pubMedExtractor.pl

           If not provided, the file 'simplePubMed.sqlite' (the
           default name used by the above script) will be used if
           found in the directory specified by -dir.

";
    warn "\n[!!] Failed to find target species:\n     $taxDat->{error}\n\n"
        if ($taxDat->{error});
    exit;
}
my ($geneMeta);

my $taxid = $taxDat->{taxid};
my $dbf   = $args->{pubmeddb};
$dbf      = "$outDir/simplePubMed.sqlite"
    if (!$dbf && -s "$outDir/simplePubMed.sqlite");


&msg("Working directory:", $tmpDir);
&msg("Output directory:", $outDir);
# &backfill_pubmed(21413195);

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

&make_metadata_file();
&ontology_pubmed();
&map_symbol();
&map_locuslink();
&map_refseq();
&ontology_go();

sub make_metadata_file {
    my $src  = "gene/DATA/gene_info.gz";
    my $vers = &_datestamp_for_file(&fetch_url($src));
    # Simple TSV file of species-specific metadata
    my $trg = sprintf("%s/Metadata-%s_Entrez-v%s.tsv", $outDir, $specID, $vers);

    unless (&output_needs_creation($trg)) {
        &msg("Using existing Metadata file:", $trg);
        return $trg;
    }
    my $tmp = "$trg.tmp";
    open (METAF, ">$tmp") || &death("Failed to write metadata file",
                                    $tmp, $!);
    
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
    print METAF join("\t", @colOut, "SymStatus") ."\n";

    my ($fh)  = &gzfh($src, $cols);
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
        my $priScore = $stat eq "O" ? "Official" : 
            $stat eq "I" ? "Interim" : "UnofficialPreferred";

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
                $priScore = "UnofficialPreferred";
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
    rename($tmp, $trg);
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

sub map_locuslink {
    my $vers = &_datestamp_for_file( &make_metadata_file() );
    my $trg = sprintf("%s/Map-%s_LocusLink-to-Entrez-v%s.mtx", 
                      $outDir, $specID, $vers);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing LocusLink map:", $trg);
        return $trg;
    }
    &capture_gene_meta();
    my @gids = sort { $a <=> $b } keys %{$geneMeta};
    my $rnum = $#gids + 1;
    my ($cnum, $nznum) = ($rnum, $rnum);
    
    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write LocusLink map", $tmp, $!);
    print MTX &_initial_mtx_block
        ( "Mapping", $rnum, $cnum, $nznum, "$specID Symbol-to-Entrez Map",
          "Simple identity matrix directly mapping LocusLink to Entrez IDs",
          "LocusLink ID", "https://www.ncbi.nlm.nih.gov/gene/?term=%s",
          "Entrez ID", $geneIdUrl, $specID);

    print MTX &_mtx_comment_block("", "This is just a very simple identity matrix that maps a LocusLink identifier to its numeric Entrez ID. Given an Entrez ID '#' the LocusLink ID would be 'LOC#'. That is, 'LOC859' is Entrez ID '859'.",
                                  "", "This mapping can of course be done very efficiently in code by simply adding or stripping 'LOC' as needed. If possible, you should take that approach. However, this matrix is generated to aid in automated analysis of diverse query inputs.", "");

    print MTX &_rowcol_meta_comment_block();

    my @meta = qw(Symbol Type Description);
    printf(MTX "%% Row %s\n", join($mtxSep, "Name", @meta));
    for my $i (0..$#gids) {
        my $id   = $gids[$i];
        my $m    = $geneMeta->{$id} || {};
        my @line = ("LOC$id", map {defined $_ ? $_ : ""} map {$m->{$_}} @meta);
        printf(MTX "%% %d %s\n", $i+1, join($mtxSep, @line));
    }

    ## Column metadata
    print MTX "% $bar\n";
    &mtx_entrez(\@gids, 'Col');

    ## The triples are just a diagonal of '1's
    print MTX &_triple_header_block( $rnum, $cnum, $nznum );
    for my $i (0..$#gids) {
        printf(MTX "%d %d 1\n", $i+1, $i+1);
    }
    close MTX;

    rename($tmp, $trg);
    &msg("Generated LocusLink to Entrez mapping", $trg);
    return $trg;
}

sub map_refseq {
    my $rsf = &make_refseq_file();
    die "WORKING";
}

sub make_refseq_file {
    my $src  = "gene/DATA/gene2refseq.gz";
    my $vers = &_datestamp_for_file(&fetch_url($src));
    my $trg = sprintf("%s/Metadata-%s_RefSeq-v%s.tsv", $outDir, $specID, $vers);
    unless (&output_needs_creation($trg)) {
        &msg("Using existing RefSeq file:", $trg);
        return $trg;
    }
    my $cols = { 
        taxid  => 'tax_id',
        gene   => 'GeneID',
        status => 'status',
        rid    => 'RNA_nucleotide_accession.version',
        pid    => 'protein_accession.version',
        mid    => 'mature_peptide_accession.version',

        # Not going to do anything with symbol here, but they're short
        # and aid in grep-debugging the file
        sym    => 'Symbol',

        # I'm not interested in GI data, they're in these columns:
        rgi    => 'RNA_nucleotide_gi',
        pgi    => 'protein_gi',
        mgi    => 'mature_peptide_gi',
    };

    my %loci;
    my ($nrow, $orow) = (0,0);

    my ($fh)  = &gzfh($src, $cols);
    my $tInd = $cols->{taxid}; # Taxid, eg 9606
    my $lInd = $cols->{gene};  # Entrez Gene ID, eg 859
    my $sInd = $cols->{status};# Status, ege REVIEWED
    my @idInds = map { $cols->{$_} } qw(rid pid mid);
    while (<$fh>) {
        s/[\n\r]+$//;
        my @row = split("\t");
        unless ($row[$tInd] == $taxid) {
            ## Stop parsing if we have already found some hits. THIS
            ## MAY BE A BAD IDEA. It should speed up parsing (neglects
            ## need to scan whole file) but there is no guarantee that
            ## taxa won't end up scattered through file in future
            ## versions:
            last if ($nrow);
            next; # Otherwise, keep scanning for taxid
        }
        map { s/^\-$// } @row; # Turn '-' cells to empty string
        my $lid  = $row[$lInd];
        my $stat = $row[$sInd];
        next unless ($lid);
        my @ids = map { $row[$_] } @idInds;
        # Remove versioning
        map { s/\.\d+$// } @ids;
        # Some species lack RNA entries (at least historically, like
        # yeast - hopefully this will eventually change). Set an empty
        # token if no RNA is seen:
        my $rid = $ids[0] || "";
        $nrow++;
        for my $p (1..2) {
            # We will collapse both "plain" and mature peptides into a
            # simple protein field.
            if (my $pid = $ids[$p]) {
                # I do not have a good sense for how the status
                # applies to RNA/Protein. Can an RNA have a different
                # status than its translated protein? We will just
                # grab the first status we encounter
                $loci{$lid}{$rid}{$pid} ||= $stat;
            }
        }
        if ($rid && !$ids[1] && !$ids[2]) {
            # RNA only, no translation
            $loci{$lid}{$rid}{""} ||= $stat;
        }
    }
    close $fh;
    my $tmp = "$trg.tmp";
    open (TSV, ">$tmp") || &death("Failed to write RefSeq file",
                                    $tmp, $!);
    my @colOut = qw(GeneID RNA Protein Status);
    print TSV join("\t", @colOut) ."\n";
    foreach my $lid (sort {$a <=> $b} keys %loci) {
        my $ldat = $loci{$lid};
        foreach my $rid (sort keys %{$ldat}) {
            my $rdat = $ldat->{$rid};
            foreach my $pid (sort keys %{$rdat}) {
                print TSV join("\t", $lid, $rid, $pid, $rdat->{$pid})."\n";
                $orow++;
            }
        }
    }
    close TSV;
    rename($tmp, $trg);
    &msg("Generated RefSeq assignments file, $orow rows from $nrow input",$trg);
    return $trg;
}


sub map_symbol {
    my $vers = &_datestamp_for_file( &make_metadata_file() );
    my $trg = sprintf("%s/Map-%s_Symbol-to-Entrez-v%s.mtx", 
                      $outDir, $specID, $vers);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing Gene Symbol map:", $trg);
        return $trg;
    }
    &capture_gene_meta();
    my (%symbols, %genes);
    my ($gnum, $nznum) = (0, 0);
    my ($lvls, $scH) = &_factorize_levels('symlevels');
    foreach my $chk (qw(Official Unofficial Unknown)) {
        die "-symlevels must include '$chk' as a level"
            unless ($scH->{uc($chk)});
    }
    my $unScore = $scH->{ uc("Unofficial") };

    # Symbols should all be the same case - but I'm not certain of
    # that. They will be collected under an upper-case key, but
    # preserve the first observed case for recording in the matrix
    foreach my $gid (sort { $a <=> $b } keys %{$geneMeta}) {
        my $meta = $geneMeta->{$gid};
        my $pri  = $meta->{Symbol} || "";
        my $ns = 0;
        my $gn;
        if ($pri) {
            # "Main" / "primary" symbol
            my $pLvl = $meta->{SymStatus} || "Unknown";
            my $pSc  = $scH->{uc($pLvl)};
            unless ($pSc) {
                $pLvl = "Unknown";
                $pSc  = $scH->{uc($pLvl)}
            }
            my $targ = $symbols{uc($pri)} ||= {
                name => $pri,
                hits => [],
            };
            $gn ||= $genes{$gid} ||= ++$gnum;
            push @{$targ->{hits}}, ($gn, $pSc);
            $nznum++;
            $ns++;
        }
        my @alis = @{$meta->{Aliases}};
        my $nUn  = $#alis + 1;
        $ns     += $nUn;

        for my $i (0..$#alis) {
            ## All aliases will have a status of Unofficial
            my $sym = $alis[$i];
            my $targ = $symbols{uc($sym)} ||= {
                name => $sym,
                hits => [],
            };
            $gn ||= $genes{$gid} ||= ++$gnum;
            push @{$targ->{hits}}, ($gn, $unScore);
            $nznum++;
        }
    }

    my @rowIds  = sort keys %symbols;
    my @gids    = sort { $genes{$a} <=> $genes{$b} } keys %genes;
    my ($rnum,$cnum) = ($#rowIds + 1, $#gids + 1);

    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write Symbol map", $tmp, $!);
    print MTX &_initial_mtx_block
        ( "Mapping", $rnum, $cnum, $nznum, "$specID Symbol-to-Entrez Map",
          "Scores are factor levels representing the nomenclature status of the assignment",
          "Gene Symbol", $symUrl,
          "Entrez ID", $geneIdUrl, $specID);


    print MTX &gene_symbol_stats_block( \%symbols, \@rowIds, $lvls, $scH );

    print MTX &_rowcol_meta_comment_block();

    ## Note row metadata
    print MTX &gene_symbol_meta( \%symbols, \@rowIds, 'Row');

    ## Column metadata
    print MTX "% $bar\n";
    print MTX &mtx_entrez(\@gids, 'Col');

    print MTX &_triple_header_block( $rnum, $cnum, $nznum );
    for my $i (0..$#rowIds) {
        my $hits = $symbols{$rowIds[$i]}{hits};
        for (my $j = 0; $j < $#{$hits}; $j += 2) {
            printf(MTX "%d %d %0.1f\n", $i+1, $hits->[$j], $hits->[$j+1]);
        }
    }
    close MTX;
    rename($tmp, $trg);
    &msg("Generated Symbol to Entrez mapping", $trg);
    return $trg;
}

sub mtx_entrez {
    my ($ids, $rc) = @_;
    &capture_gene_meta();
    return &_generic_meta_block($id, $rc, $geneMeta, 
                                [qw(Symbol Type Description)]);
}

sub ontology_pubmed {
    my $src  = "gene/DATA/gene2pubmed.gz";
    my $vers = &_datestamp_for_file(&fetch_url($src));
    my $trg = sprintf("%s/Ontology-%s_Entrez-to-PubMed-v%s.mtx", 
                      $outDir, $specID, $vers);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing PubMed file:", $trg);
        return $trg;
    }
    unless ($dbf) {
        &err("Can not create PubMed ontology",
             "Please provide SQLite DB path with -pubmeddb",
             "The DB can be generated with pubMedExtractor.pl");
        return "";
    }

    &capture_gene_meta();
    &msg("Parsing PubMed");

    ## If you see 'expected column' errors, you will need to change
    ## the right hand value (after the '=>') of the offending column,
    ## after determining what the new column name is:
    my $cols = { 
        taxid  => 'tax_id',
        gene   => 'GeneID',
        pmid   => 'PubMed_ID',
    };
    my ($fh) = &gzfh($src, $cols);
    my $tInd = $cols->{taxid}; # Taxid, eg 9606
    my $oInd = $cols->{pmid};  # PubMed ID, eg 9873079
    my $lInd = $cols->{gene};  # Entrez Gene ID, eg 859
    my ($nznum, $nont) = (0,0);
    my (%genes, %ontMeta);

    while (<$fh>) {
        s/[\n\r]+$//;
        my @row = split("\t");
        unless ($row[$tInd] == $taxid) {
            ## Stop parsing if we have already found some hits. THIS
            ## MAY BE A BAD IDEA. It should speed up parsing (neglects
            ## need to scan whole file) but there is no guarantee that
            ## taxa won't end up scattered through file in future
            ## versions:
            last if ($nont);
            next; # Otherwise, keep scanning for taxid
        }
        map { s/^\-$// } @row; # Turn '-' cells to empty string
        my $oid = $row[$oInd];
        my $om = $ontMeta{$oid} ||= {
            order       => ++$nont,
        };
        # Capture orderID/score pairs
        push @{$genes{ $row[ $lInd ] }}, ( $om->{order}, 1 );
        $nznum++;
    }
    close $fh;
    warn "     ... writing matrix market file ...\n";

    my @ontIds = sort { $ontMeta{$a}{order} <=> 
                            $ontMeta{$b}{order} } keys %ontMeta;
    my @gids   = sort { $a <=> $b } keys %genes;
    my ($rnum, $cnum) = ($#gids + 1, $#ontIds + 1);

    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write PubMed ontology", $tmp, $!);
    print MTX &_initial_mtx_block
        ( "Ontology", $rnum, $cnum, $nznum, "$specID Entrez PubMed",
          "PubMed assignments for $specID Entrez genes",
          "Entrez ID", $geneIdUrl,
          "PubMed ID", "https://www.ncbi.nlm.nih.gov/pubmed/%s", $specID);
    print MTX "%
% These data were extracted from the NLM/NCBI FTP sites:
%   https://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz
%   https://ftp.ncbi.nih.gov/pubmed/baseline/
%   https://ftp.ncbi.nih.gov/pubmed/updatefiles/
";
    print MTX &_setfisher_comment_block();
    print MTX &_min_set_mtx_block( 5 );
    print MTX &_max_set_mtx_block
        ( 15, "An exhaustive survey of the Drop Bear transcriptome" );
    print MTX &_min_onto_mtx_block( 2 );
    print MTX &_rowcol_meta_comment_block();

    ## Row metadata
    &mtx_entrez(\@gids, 'Row');
    print MTX "% $bar\n";

    ## Column metadata
    
    &msg("Extracting PubMed metadata", $dbf);
    &_sqlite_tools();

    my @meta = qw(Date Description);
    printf(MTX "%% Col %s\n", join($mtxSep, "Name", @meta));
    for my $i (0..$#ontIds) {
        my $oId = $ontIds[$i];
        my ($vers, $dt, $desc) = &_pmid_from_db( $oId );
        ## Not all PMIDs are in the XML files (eg Book entries)
        ## Use NCBI EFetch to recover these:
        ($vers, $dt, $desc) = &backfill_pubmed( $oId) unless ($desc);
        printf(MTX "%% %d %s\n", $i+1, join($mtxSep, $oId, $dt, $desc));
    }

    &msg("Writing connections");
    print MTX &_triple_header_block( $rnum, $cnum, $nznum );
    for my $i (0..$#gids) {
        my $hits = $genes{$gids[$i]};
        for (my $j = 0; $j < $#{$hits}; $j += 2) {
            printf(MTX "%d %d %d\n", $i+1, $hits->[$j], $hits->[$j+1]);
        }
    }
    close MTX;
    rename($tmp, $trg);
    &msg("Generated Entrez PubMed ontology", $trg);
    return $trg;
}

sub ontology_go {
    my $src  = "gene/DATA/gene2go.gz";
    my $vers = &_datestamp_for_file(&fetch_url($src));
    my $trg = sprintf("%s/Ontology-%s_Entrez-to-GeneOntology-v%s.mtx", 
                      $outDir, $specID, $vers);
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

    my ($lvls, $scH) = &_factorize_levels('eclevels', 'upper');

    my ($fh) = &gzfh($src, $cols);
    my $tInd = $cols->{taxid}; # Taxid, eg 9606
    my $gInd = $cols->{goid};  # GO ID, eg GO:0005886
    my $eInd = $cols->{ec};    # Evidence code, eg TAS
    my $lInd = $cols->{gene};  # Entrez Gene ID, eg 859
    my (%genes, %ontMeta);
    my ($nznum, $nont) = (0,0);

    while (<$fh>) {
        s/[\n\r]+$//;
        my @row = split("\t");
        unless ($row[$tInd] == $taxid) {
            ## Stop parsing if we have already found some hits. THIS
            ## MAY BE A BAD IDEA. It should speed up parsing (neglects
            ## need to scan whole file) but there is no guarantee that
            ## taxa won't end up scattered through file in future
            ## versions:
            last if ($nont);
            next; # Otherwise, keep scanning for taxid
        }
        map { s/^\-$// } @row; # Turn '-' cells to empty string
        my $gid = $row[$gInd];
        ## Capture term metadata once
        my $gm = $ontMeta{$gid} ||= {
            Description => $row[ $cols->{term} ],
            Category    => $row[ $cols->{cat}  ],
            order       => ++$nont,
        };
        my $sc = $scH->{ uc($row[ $eInd ]) };
        unless ($sc) {
            unless (defined $sc) {
                &err("Unknown Evidence code '$row[ $eInd ]'",
                     "This code was not in -eclevels, it will be ignored");
                $scH->{ uc($row[ $eInd ]) } = 0;
            }
            next;
        }
        # Capture orderID/score pairs
        push @{$genes{ $row[ $lInd ] }}, ( $gm->{order}, $sc );
        $nznum++;
    }
    close $fh;
    warn "     ... writing matrix market file ...\n";

    my @ontIds = sort { $ontMeta{$a}{order} 
                        <=> $ontMeta{$b}{order} } keys %ontMeta;
    my @gids   = sort { $a <=> $b } keys %genes;
    my ($rnum, $cnum) = ($#gids + 1, $#ontIds + 1);

    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write GO ontology", $tmp, $!);
    print MTX &_initial_mtx_block
        ( "Ontology", $rnum, $cnum, $nznum, "$specID Entrez GeneOntology",
          "GeneOntology assignments for $specID Entrez genes",
          "Entrez ID", $geneIdUrl,
          "GO Term", "http://amigo.geneontology.org/amigo/term/%s", $specID);
    print MTX &_setfisher_comment_block();
    print MTX &_min_set_mtx_block( 7 );
    print MTX &_max_set_mtx_block( 5, "Catalytic process" );
    print MTX &_min_onto_mtx_block( 2 );
    print MTX &_factor_map_block
        ( $lvls, $scH, "Evidence Codes",
          ["Lower factor values should represent lower confidence evidence",
           "http://geneontology.org/page/guide-go-evidence-codes"]);
    print MTX &_rowcol_meta_comment_block();

    ## Row metadata
    &mtx_entrez(\@gids, 'Row');
    print MTX "% $bar\n";

    ## Column metadata
    my @meta = qw(Description Category);
    printf(MTX "%% Col %s\n", join($mtxSep, "Name", @meta));
    for my $i (0..$#ontIds) {
        my $oId = $ontIds[$i];
        my $dat = $ontMeta{$oId};
        printf(MTX "%% %d %s\n", $i+1, join($mtxSep, $oId,
                                            map { $dat->{$_} } @meta));
    }

    print MTX &_triple_header_block( $rnum, $cnum, $nznum );

    for my $i (0..$#gids) {
        my $hits = $genes{$gids[$i]};
        for (my $j = 0; $j < $#{$hits}; $j += 2) {
            printf(MTX "%d %d %d\n", $i+1, $hits->[$j], $hits->[$j+1]);
        }
    }
    close MTX;
    rename($tmp, $trg);
    &msg("Generated Entrez GO ontology", $trg);
    return $trg;
}

sub extract_taxa_info {
    ## Used to deconvolute user species request into a formal
    ## species. In particular, we'll need the taxid (eg 9606 for
    ## human) to extract relevant subsets of information.
    my $req = shift;
    return { error => "No taxa request specified" }  unless ($req);
    my $srcFile = "$tmpDir/names.dmp";
    if (&source_needs_recovery($srcFile)) {
        my $tgz = &fetch_url("pub/taxonomy/taxdump.tar.gz");
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

sub _setfisher_comment_block {
    return "% $bar
% The filters described below can be recognized by the SetFisher
% package as pruning parameters. The filters are applied recursively,
% removing terms, then genes, and then repeating until no further
% alterations occur to the matrix.
% $bar
";
}

sub _min_onto_mtx_block {
    my $def = shift;
    $def = 2 unless (defined $def);
    return "%
% A major distortion in Fisher-based testing can come from
% 'questionable' genes. These include real-but-untranscribed entities
% (eg pseudogenes), speculative entries that will eventually be
% removed from the transcriptome, or rare genes that are expressed
% only in certain tissues or developmental stages. From a
% marbles-in-an-urn perspective, all these categories represent
% marbles that can NEVER be removed from the urn. The effect is
% generally to inflate the significance of those that you do, since
% the world size is larger than it really should be AND these
% unselectable marbles generally have sparse annotation compared to
% the 'real' ones that can be selected. This reduced level of
% annotation is used by minOntoSize to remove genes that are likely to
% be speculative. The filter below removes genes with fewer than the
% indicated number of terms assigned to it.
%
%% DEFAULT minOntoSize $def
";
}

sub _min_set_mtx_block {
    my $def = shift;
    $def = 5 unless (defined $def);
    return "%
% Terms with few genes assigned to them struggle to reach statistical
% significance. Excluding them removes some distractions, but also
% helps minimize multiple-testing penalties on 'long shot' terms. The
% threshold below represents the minimum number of genes to be
% assigned to a term for the term to be kept.
%
%% DEFAULT minSetSize $def
";

}

sub _max_set_mtx_block {
    my ($def, $example) = @_;
    $def = 10 unless (defined $def);
    return "%
% Terms with many genes tend to be less useful in biological
% interpretation - for example:
%  '$example'
% We have also observed that they tend to be disproportionately
% significant as they are greatly impacted by the effect mentioned
% under minOntoSize (un-selectable genes distorting the significance
% of enriched sets). If a term is assigned to a greater percentage of
% the world than the value listed below it will be excluded from the
% matrix used in analysis. In addition to removing clutter from
% reports, it brings a small multiple testing benefit.
%
%% DEFAULT maxSetPerc $def
";
}

sub backfill_pubmed {
    my $pmid = shift;
    # my $url  = "https://www.ncbi.nlm.nih.gov/pubmed/?report=xml&format=json&term=$pmid";
    # my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ecitmatch.cgi?db=pubmed&retmode=json&id=$pmid";
    my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id=$pmid";
    my $file = &get_url($url, "BackFill-$pmid.xml");
    unless (-s $file) {
        &err("Failed to recover XML data for PubMed", $url);
        return ( -2, "", "-Error recovering XML from PubMed-" );
    }

    my $twig = XML::Twig->new();
    $twig->parsefile( $file );
    my $doc  = $twig->root();
    unless ($doc) {
        &err("PubMed XML document is ... odd", $file);
        return ( -2, "", "-Error parsing XML from PubMed-" );
    }
    $doc = $doc->first_child();
    $doc = $doc->first_child() if ($doc);
    unless ($doc) {
        &err("PubMed does not appear to have information on PMID:$pmid", $file);
        my $dat = [-2, "", "-This article is apparently no longer at PubMed-" ];
        &msg("PubMed ID $pmid appears to have no information in PubMed");
        &_backfill_tools();
        $bfClear->execute( $pmid );
        $bfSet->execute( $pmid, $dat->[1], $dat->[2]);
        return @{$dat};
    }
    my @chk  = $doc->get_xpath('PMID');
    my $dat  = [0,'',''];
    if ($#chk == 0 && $chk[0]->text() == $pmid) {
        ## Double-check that we got the right XML snippet
        my @dts = $doc->get_xpath('//PubDate');
        if ($#dts > -1) {
            # Found a date block
            my $y = [ map { $_->text() } $dts[0]->get_xpath('Year') ];
            my $m = [ map { $_->text() } $dts[0]->get_xpath('Month') ];
            my $d = [ map { $_->text() } $dts[0]->get_xpath('Day') ];
            $dat->[1] = &parse_pubmed_date($pmid, $y, $m, $d);
        }
        foreach my $di (qw(ArticleTitle VernacularTitle BookTitle Abstract)) {
            my @nodes = $doc->get_xpath("//$di");
            if (my $found = $nodes[0]) {
                if (my $txt = $found->text()) {
                     $txt = substr($txt, 0, $maxAbst).'...' 
                         if (length($txt) > $maxAbst && $di eq 'Abstract');
                     $dat->[2] = $txt;
                     last;
                }
            }
        }
        ## Add the data to the database
        unless ($dat->[2]) {
            $dat->[2] = "-Title not identified in XML-";
            &msg("Failed to surgically recover title for PMID $pmid");
        }
        &_backfill_tools();
        $bfClear->execute( $pmid );
        $bfSet->execute( $pmid, $dat->[1], $dat->[2]);
    } else {
        $dat = [-1, '', 'Error recovering XML information from PubMed'];
    }
    return @{$dat};
}

sub _sqlite_tools {
    return $dbh if ($dbh);
    ## Read-only connection: https://stackoverflow.com/a/34360995
    $dbh = DBI->connect("dbi:SQLite:dbname=$dbf",'','', {
        sqlite_open_flags => DBD::SQLite::OPEN_READONLY,
        AutoCommit => 1,
        RaiseError => 1,
        PrintError => 0 });
    $getPMID = $dbh->prepare
        ("SELECT vers, pubdate, title FROM pmid WHERE pmid = ?");
    return $dbh;
}

sub _pmid_from_db {
    my $oId = shift || 0;
    $getPMID->execute( $oId );
    my $dat = $getPMID->fetchall_arrayref();
    my ($newest) = sort { $b->[0] <=> $a->[0] } @{$dat};
    return @{$newest || [0, '','']};
}

sub _backfill_tools {
    &_sqlite_tools( @_ );
    return $bfdbh if ($bfdbh);
    $bfdbh = DBI->connect("dbi:SQLite:dbname=$dbf",'','', {
        AutoCommit => 1,
        RaiseError => 1,
        PrintError => 0 });
    $bfClear = $bfdbh->prepare("DELETE FROM pmid WHERE pmid = ? AND vers = 0");
    $bfSet   = $bfdbh->prepare
        ("INSERT INTO pmid (pmid, vers, xsid, pubdate, title)".
         " VALUES (?, 0, 0, ?, ?)");
    return $bfdbh;
}
