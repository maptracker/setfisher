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
        Type        => 'type_of_gene',
        Status      => 'Nomenclature_status',
        Aliases     => 'Synonyms',
        Description => 'description',
    };
    my @colOut = qw(GeneID Symbol Type Status Description Aliases);
    print METAF join("\t", @colOut) ."\n";

    my ($fh)  = &gzfh("gene/DATA/gene_info.gz", "gene_info.gz", $cols);
    my @getInds = map { $cols->{$_} } @colOut;
    my $tInd = $cols->{TaxID};   # Taxid, eg 9606
    my $lInd = $cols->{GeneID};  # GeneID, eg 859
    my $aInd = $cols->{Aliases}; # Alt symbols, eg LGMD1C|LQT9|VIP-21|VIP21
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
        # I am not certain if it is needed, but I want to clean up the
        # alias field in case there is odd whitespace or null entries.
        my $aliases = $row[ $aInd ] || "";
        $aliases =~ s/^\s+//; $aliases =~ s/\s+$//; # Edge whitespace
        $aliases =~ s/\|[\s\|]+/\|/g; # Null entries
        $aliases = "" if ($aliases eq '|');
        $row[ $aInd ] = $aliases;
        ## Add the entry to the TSV file:
        print METAF join("\t", map { $row[$_] } @getInds)."\n";
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
    $head =~ s/\[\n\r]+$//;
    my @cols = split(/\t/, $head);
    while (<METAF>) {
        s/\[\n\r]+$//;
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
    my @gids = sort { $a <=> $b } keys %{$geneMeta};
    my (%symbols, %genes);
    my $gnum = 0;
    my %symSc = ( O => 1,   I => 0.5, "" => 0.3 );
    foreach my $gid (sort { $a <=> $b } keys %{$geneMeta}) {
        my $meta = $geneMeta->{$gid};
        my $stat = $meta->{Status} || "";

    }

    die "Working"
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
    my ($nrow, $ngo) = (0,0);

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
        $nrow++;
    }
    close $fh;
    warn "     ... writing matrix market file ...\n";

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
