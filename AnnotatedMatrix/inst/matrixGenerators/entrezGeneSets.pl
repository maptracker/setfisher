#!/usr/bin/perl -w

use strict;
my $scriptDir;
our $defTmp = "/tmp/entrezGeneSets";
my $defFtp  = "ftp.ncbi.nih.gov";
my $defSym  = "Unknown,Unofficial,UnofficialPreferred,Interim,Official";
my $defEC   = "P,IEA,NAS,IRD,IBD,IBA,RCA,IGC,ISS,ISA,ISM,ISO,TAS,EXP,IEP,IPI,IMP,IGI,IDA";

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
our ($args, $outDir, $clobber, $tmpDir, $maxAbst, $bar);
my ($dbh, $getPMID, $bfdbh, $bfClear, $bfSet, %tasks);

use DBI;
use DBD::SQLite;
use Archive::Tar;
use XML::Twig;

my $sources = {
    GeneInfo    => "gene/DATA/gene_info.gz",
    Gene2RefSeq => "gene/DATA/gene2refseq.gz",
    Gene2PubMed => "gene/DATA/gene2pubmed.gz",
    Gene2GO     => "gene/DATA/gene2go.gz",
    Taxonomy    => "pub/taxonomy/taxdump.tar.gz",
    GO          => "ftp://ftp.geneontology.org/go/ontology/go.obo",
};

my $vers       = &fetch_all_files();
my $taxDat     = &extract_taxa_info( $args->{species} || $args->{taxa} );

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
my ($geneMeta, $goMeta);
my @stndMeta  = qw(Symbol Type Description);
my $nsUrl = {
    Symbol      => 'https://www.ncbi.nlm.nih.gov/gene/?term=%s%%5Bsym%%5D',
    EntrezGene  => 'https://www.ncbi.nlm.nih.gov/gene/%s', # Integer IDs
    LocusLink   => 'https://www.ncbi.nlm.nih.gov/gene/?term=%s', # For LOC###
    PubMed      => 'https://www.ncbi.nlm.nih.gov/pubmed/%s',
};

my $taxid = $taxDat->{taxid};
my $dbf   = $args->{pubmeddb};
$dbf      = "$outDir/simplePubMed.sqlite"
    if (!$dbf && -s "$outDir/simplePubMed.sqlite");

my $auth     = "Entrez";
my $authLong = "$auth ## Data repository at the National Center for Biotechnology Information";



&msg("'Release' version: $vers");
&msg("Working directory:", $tmpDir);
&msg("Output directory:",  $outDir);
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



## Stash deduplicated file store - not on all systems
my $stashMeta = {
    Authority  => $auth,
    Version    => $vers,
    MatrixType => "Map",
    FileType   => "AnnotatedMatrix",
    Format     => "MatrixMarket",
};




&make_entrez_metadata_file();
&map_symbol();
&map_locuslink();
&ontology_go();
die;
&map_refseq();
&ontology_pubmed();

sub fetch_all_files {
    ## The files being used here are not part of a unified
    ## release. Each specific source file will be versioned with its
    ## date stamp from the FTP server that it was recovered
    ## from. However, we are also organizing the files under a single
    ## version. For this reason, today's date will be taken as the
    ## "release version", and all files will be bulk recovered at once
    ## to try to get a semblance of consistency.
    &msg("Downloading source files");
    foreach my $tag (sort keys %{$sources}) {
        my $src   = $sources->{$tag};
        my $fVers = &_datestamp_for_file(&fetch_url($src));
        &msg("    $tag :  $fVers");
    }
    my ($sec,$min,$hour,$mday,$mon,$year) = localtime();
    my $vFile = "$tmpDir/BuildDate.txt";
    unless (-s $vFile && !$clobber) {
        open(VF, ">$vFile") || &death
            ("Failed to write date to version file", $vFile, $!);
        printf(VF "%04d-%02d-%02d\n", $year+1900, $mon+1, $mday);
        close VF;
    }
    open(VF, "<$vFile") || &death
        ("Failed to read date from version file", $vFile, $!);
    my $rv = <VF>;
    $rv =~ s/[\n\r]+$//;
    close VF;
    return $rv;
}

sub make_entrez_metadata_file {
    my $src   = $sources->{GeneInfo};
    my $fVers = &_datestamp_for_file(&fetch_url($src));
    # Simple TSV file of species-specific metadata
    my $trg = sprintf("%s/Metadata@%s-%s@%s@%s.tsv",
                      &primary_folder($auth, $vers), $specID,
                      "EntrezGene", $auth, $fVers);

    unless (&output_needs_creation($trg)) {
        &msg("Using existing Metadata file:", $trg) unless ($tasks{MakeMeta}++);
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

sub make_go_metadata_file {
    my $src   = $sources->{GO};
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

    # Simple TSV file of GeneOntology metadata
    my $trg = sprintf("%s/Metadata@%s@%s.tsv",
                      &primary_folder($auth, $vers), "GeneOntology", $fVers);

    ## http://www.geneontology.org/page/ontology-relations
    my $transitive = {
        "is_a" => {
            "is_a"      => "is_a",
            "part_of"   => "part_of",
            "regulates" => "regulates",
        },
        "part_of" => {
            "part_of"   => "part_of",
            "is_ia"     => "part_of",
        },
        "regulates" => {
            "is_a"      => "regulates",
            "part_of"   => "regulates",
        },
    };
    ## ignoring has_part : inverse of part_of
    ## ignoring negatively_regulates and positively_regulates

    unless (&output_needs_creation($trg)) {
        &msg("Using existing Metadata file:", $trg) unless ($tasks{GoMeta}++);
        return $trg;
    }

    &msg("Parsing basic GeneOntology metadata");
    my $tmp = "$trg.tmp";
    open (METAF, ">$tmp") || &death("Failed to write metadata file", $tmp, $!);
    my @head = qw(id name is_obsolete namespace is_a part_of intersection_of regulates synonym def);
    my $hMap = {
        id => "ID",
        name => "Name",
        namespace => "Namespace",
        def       => "Definition",
        synonym   => "Synonyms",
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
                    ## "downsampling" positive and negative regulation
                    $k = "regulates" if ($k =~ /_regulates/);
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
    &msg("Generated Entrez metadata file", $trg);
    return $trg;
}

sub _go_meta_record {
    my ($rec, $fh, $head, $nextK) = @_;
    my $k = $rec->{now};
    my $d = $rec->{dat};
    if ($k eq 'Term') {
        if (my $id = $d->{id}) {
            my @row;
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
    die "WORKING HERE";
}

sub capture_gene_meta {
    return $geneMeta if ($geneMeta);
    $geneMeta  = {};
    my $mdFile = &make_entrez_metadata_file();
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
    my $fVers = &_datestamp_for_file( &make_entrez_metadata_file() );
    my ($nsi, $nsj) = ("LocusLink", "EntrezGene");
    my @fbits = ("Map", $specID, $nsi, $nsj, $auth, $fVers);
    my $trg   = &primary_path(@fbits);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing $nsi map:", $trg);
        return $trg;
    }
    &capture_gene_meta();
    my @gids = sort { $a <=> $b } keys %{$geneMeta};
    my $rnum = $#gids + 1;
    my ($cnum, $nznum) = ($rnum, $rnum);
    
    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write $nsi map", $tmp, $!);
    print MTX &_initial_mtx_block
        ( "Mapping", $rnum, $cnum, $nznum, "$auth $specID $nsi-to-$nsj Map",
          "Simple identity matrix directly mapping $specID $nsi to $nsj IDs",
          "All scores are 1",
          $nsi, $nsj);

    print MTX &_dim_block({
        RowDim    => $nsi,
        RowUrl    => $nsUrl->{$nsi},
        ColDim    => $nsj,
        ColUrl    => $nsUrl->{$nsj}, 
        Authority => $authLong,
        Source    => &ftp_url($sources->{GeneInfo}),
        Version   => $vers });


    print MTX &_citation_MTX();
    print MTX &_species_MTX( $taxDat->{'scientific name'} );

    print MTX &_mtx_comment_block("", "This is just a very simple identity matrix that maps a LocusLink identifier to its numeric Entrez ID. Given an Entrez ID '#' the LocusLink ID would be 'LOC#'. That is, 'LOC859' is Entrez ID '859'.",
                                  "", "This mapping can of course be done very efficiently in code by simply adding or stripping 'LOC' as needed. If possible, you should take that approach. However, this matrix is generated to aid in automated analysis of diverse query inputs.", "");

    print MTX &_rowcol_meta_comment_block();

    ## Note row metadata
    ## Locus IDs are just the GeneIDs with "LOC" prefixed.
    my @lids = map { "LOC$_"} @gids;
    ## Copy over metadata:
    map { $geneMeta->{"LOC$_"} = $geneMeta->{$_} } @gids;

    print MTX &_generic_meta_block(\@lids, 'Row', $geneMeta, \@stndMeta);

    ## Column metadata
    print MTX "% $bar\n";
    print MTX &_generic_meta_block(\@gids, 'Col', $geneMeta, \@stndMeta);

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
    my $src   = $sources->{Gene2RefSeq};
    my $fVers = &_datestamp_for_file(&fetch_url($src));
    my $trg = sprintf("%s/Metadata-%s_RefSeq-v%s.tsv", $outDir, $specID, $fVers);
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
    my $fVers = &_datestamp_for_file( &make_entrez_metadata_file() );
    my ($nsi, $nsj) = ("Symbol", "EntrezGene");
    my @fbits = ("Map", $specID, $nsi, $nsj, $auth, $fVers);
    my $trg   = &primary_path(@fbits);
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
        ( "Mapping", $rnum, $cnum, $nznum, "$auth $specID $nsi-to-$nsj Map",
          "Accession conversion from $nsi to $nsj",
          "Factor levels representing nomenclature status of the assignment",
          $nsi, $nsj);

    print MTX &_dim_block({
        RowDim    => $nsi,
        RowUrl    => $nsUrl->{$nsi},
        ColDim    => $nsj,
        ColUrl    => $nsUrl->{$nsj}, 
        Authority => $authLong,
        Source    => join('/', "ftp:/", $args->{ftp}, $sources->{GeneInfo}),
        Version   => $vers });

    print MTX &_citation_MTX();
    print MTX &_species_MTX( $taxDat->{'scientific name'} );
    print MTX &_filter_block({
        TossLevel => "[,][Unknown] ## 'Unknown' entries are deprecated assignments."
                             });
    print MTX &gene_symbol_stats_block( \%symbols, \@rowIds, $lvls, $scH );

    print MTX &_rowcol_meta_comment_block();

    ## Note row metadata
    print MTX &gene_symbol_meta( \%symbols, \@rowIds, 'Row');

    ## Column metadata
    print MTX "% $bar\n";
    print MTX &_generic_meta_block(\@gids, 'Col', $geneMeta, \@stndMeta);

    print MTX &_triple_header_block( $rnum, $cnum, $nznum );
    for my $i (0..$#rowIds) {
        my $hits = $symbols{$rowIds[$i]}{hits};
        for (my $j = 0; $j < $#{$hits}; $j += 2) {
            printf(MTX "%d %d %0.1f\n", $i+1, $hits->[$j], $hits->[$j+1]);
        }
    }
    close MTX;
    rename($tmp, $trg);
    &msg("Generated $nsi to $nsj mapping", $trg);
    return $trg;
}

sub ontology_pubmed {
    my $src   = $sources->{Gene2PubMed};
    my $fVers = &_datestamp_for_file(&fetch_url($src));
    my ($nsi, $nsj) = ("EntrezGene", "PubMed");
    my @fbits = ("Ontology", $specID, $nsi, $nsj, $auth, $fVers);
    my $trg   = &primary_path(@fbits);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing $nsj file:", $trg);
        return $trg;
    }
    unless ($dbf) {
        &err("Can not create $nsj ontology",
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
        ( "Ontology", $rnum, $cnum, $nznum, "$auth $specID $nsi $nsj Ontology",
          "$nsj assignments for $specID $nsi",
          "Scores are all simply 1",
          $nsi, $nsj);

    print MTX &_dim_block({
        RowDim    => $nsi,
        RowUrl    => $nsUrl->{$nsi},
        ColDim    => $nsj,
        ColUrl    => $nsUrl->{$nsj}, 
        Authority => $authLong,
        Source    => ["https://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz",
                      "https://ftp.ncbi.nih.gov/pubmed/baseline/",
                      "https://ftp.ncbi.nih.gov/pubmed/updatefiles/"],
        Version   => $vers });

    print MTX &_citation_MTX();
    print MTX &_species_MTX( $taxDat->{'scientific name'} );

    print MTX &_filter_block({
        MinColCount => "5 ## Articles with few genes assigned to them struggle to reach statistical significance",
        MaxColCount => "15% ## Articles containing a large fraction of the genome are rarely informative",
        MinRowCount => "2 ## Genes with few publications generally do not bring much insight to the analysis",
                             });

    print MTX &_rowcol_meta_comment_block();

    ## Row metadata
    print MTX &_generic_meta_block(\@gids, 'Row', $geneMeta, \@stndMeta);

    print MTX "% $bar\n";

    ## Column metadata
    
    &msg("Extracting PubMed metadata", $dbf);
    &_sqlite_tools();
    my $pmMeta = {};
    foreach my $oId (@ontIds) {
        my ($pmVers, $dt, $desc) = &_pmid_from_db( $oId );
        ## Not all PMIDs are in the XML files (eg Book entries)
        ## Use NCBI EFetch to recover these:
        ($pmVers, $dt, $desc) = &backfill_pubmed( $oId) unless ($desc);
        $pmMeta->{$oId} = {
            Date        => $dt,
            Description => $desc,
        };
    }
    print MTX &_generic_meta_block(\@ontIds, 'Col', $pmMeta,
                                   [qw(Date Description)]);

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
    my $src   = $sources->{Gene2GO};
    my $fVers = &_datestamp_for_file(&fetch_url($src));
    my ($nsi, $nsj) = ("EntrezGene", "GeneOntology");
    my @fbits = ("Ontology", $specID, $nsi, $nsj, $auth, $fVers);
    my $trg   = &primary_path(@fbits);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing $nsj file:", $trg);
        return $trg;
    }
    &capture_gene_meta();
    &go_meta();
    die "WORK HERE TOO";
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
        ( "Ontology", $rnum, $cnum, $nznum, "$auth $specID $nsi-to-$nsj Ontology",
          "$specID $nsi $nsj ontology, as assigned by $auth",
          "Scores are factors representing the evidence code for the assignment",
          $nsi, $nsj);

    print MTX &_dim_block({
        RowDim    => $nsi,
        RowUrl    => $nsUrl->{$nsi},
        ColDim    => $nsj,
        ColUrl    => $nsUrl->{$nsj}, 
        Authority => $authLong,
        Source    => &ftp_url($src),
        Version   => $vers });

    print MTX &_citation_MTX();
    print MTX &_species_MTX( $taxDat->{'scientific name'} );

    print MTX &_filter_block({
        MinColCount => "7 ## Terms with few genes assigned to them struggle to reach statistical significance",
        MaxColCount => "10% ## Terms covering a large fraction of the genome are rarely informative. Note this will exclude many parent terms",
        MinRowCount => "2 ## Genes with few terms assigned to them will not bring much insight to the analysis",
        TossLevel => "[,][P,IEA] ## P is an ancient evidence code, IEA tend to have high false positives (but may represent many of your assignments)"
                             });

    print MTX &_factor_map_block
        ( $lvls, $scH, "Evidence Codes",
          ["Lower factor values generally represent lower confidence evidence",
           "http://geneontology.org/page/guide-go-evidence-codes"]);

    print MTX &_rowcol_meta_comment_block();

    ## Row metadata
    print MTX &_generic_meta_block(\@gids, 'Row', $geneMeta, \@stndMeta);
    print MTX "% $bar\n";

    ## Column metadata
    print MTX &_generic_meta_block(\@ontIds, 'Col', \%ontMeta, 
                                   [qw(Description Category)]);

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
    return { error => "-species is not specified" }  unless ($req);
    my $srcFile = "$tmpDir/names.dmp";
    if (&source_needs_recovery($srcFile)) {
        my $tgz = &fetch_url($sources->{Taxonomy});
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
    return { error => "-species '$req' could not be found in $srcFile" };
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

sub _citation_MTX {
    return "%\n". &_default_parameter( "Citation", "Gene [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; 2004 – [$vers]. Available from: https://www.ncbi.nlm.nih.gov/gene/")."%\n";
}
