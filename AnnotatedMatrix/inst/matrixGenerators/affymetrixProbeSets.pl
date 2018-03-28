#!/usr/bin/perl -w

use strict;
my $scriptDir;
our $defTmp = "/tmp/AffymetrixGeneSets";
my $defFtp  = "www.affymetrix.com";
my $defCtx  = "unknown,downstream,upstream,intron,splice-site,exon,UTR-3,UTR-5,cds,synon,missense,nonsense";


our $defaultArgs = {
    ftp            => $defFtp,
    contextlevels  => $defCtx,
    dir            => ".",
    array          => "UNKNOWN",
    version        => "",
    bycontext      => "",
};

BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}

use lib "$scriptDir";
require Utils;
use IO::Uncompress::Unzip qw(unzip);
use Text::CSV;

our ($args, $outDir, $clobber, $tmpDir, $maxAbst, $bar);
my ($dbh, $getPMID, $bfdbh, $bfClear, $bfSet, %tasks);

my @stndMeta   = qw(Symbol Description);
my %defColDef = ( Symbol => "Entrez Gene symbol as provided by Affymetrix",
                  Description => "Short descriptive text for an Entrez Gene" );

my $nsUrl = {
    EntrezGene     => 'https://www.ncbi.nlm.nih.gov/gene/%s', # Integer IDs
};


my $aReq = $args->{array} || "";
my $vers = $args->{version} || "";

if ($vers =~ /(\d+)/) { $vers = $1; }

my ($array, $arrayToken);

if ($aReq =~ /(gw.?6|genomewide.*6)/i) {
    # http://www.affymetrix.com/Auth/analysis/downloads/na34/genotyping/GenomeWideSNP_6.cn.na34.annot.csv.zip
    $array      = "GenomeWideSNP_6";
    $arrayToken = "GW6";
    $vers ||= 34; # Most recent (Sep 2013) as of March 2018
} elsif ($aReq) {
    warn "Array request '$aReq' not recognized\n";
}

$nsUrl->{AffyProbeSet} = "https://www.affymetrix.com/analysis/netaffx/mappingfullrecord.affx?pk=${array}:%s";

my $naRegister = "http://www.affymetrix.com/site/terms.affx?dest=register&buttons=on";

if (!$array || !$vers || $args->{h} || $args->{help}) {
    warn "
Usage:

This program will generate gene sets for the SetFisher enrichment
package. It will recover information from the Affymetrix web site
based on an array design and a NetAffx version that you provide.

Required Argument:

    -array The Affymetrix array design desired

Optional Arguments:

   -tmpdir Default '$defTmp'. Directory holding downloaded files

      -dir Default '.'.  Output directory that will contain generated
           matrix and metadata files

  -clobber Default 0, which will preserve any already-generated
           files. Set to 1 to regenerate matrix files, and any value
           greater than 1 to also re-download base files from Entrez.

     -name A name to use for the files. By default this will be chosen
           from available taxonomy information, preferring the Genbank
           Common Name

 -contextlevels Default '$defCtx'. This series maps the 'context' of
           the probe/marker relative to the nearby gene. Used by
           genomic arrays.

  -bycontext Default '' (empty string), which will cause genomic
           matrices to be scored by the distance from the probe to the
           gene. If a non-false value, then the matrix will instead be
           scored as a factor based on the context.

Please also review the Affymetrix terms of use here:

  $naRegister

You will likely need to download the files manually

";
    exit;
}

my $auth     = "Affymetrix";
my $authLong = "$auth ## Purveyor of nucleotide hybridization profiling systems. http://www.affymetrix.com/";
my $versToken = "NA$vers";

## Stash deduplicated file store - not on all systems
my $stashMeta = {
    Authority  => $auth,
    'Version'  => $vers,
    MatrixType => "Map",
    FileType   => "AnnotatedMatrix",
    Format     => "MatrixMarket",
};


&msg("'Release' version: $vers");
&msg("Output directory:",  $outDir);


my $baseUrl = "http://www.affymetrix.com/Auth/analysis/downloads/na" . $vers;


if ($arrayToken eq 'GW6') {
    &parse_GenomeWide();
} else {
    die "Have not built URL logic for non-GW6 arrays yet"
}

sub parse_GenomeWide {
    ## Looking for files like:
    ##    GenomeWideSNP_6.na34.annot.csv.zip
    ##    GenomeWideSNP_6.cn.na34.annot.csv.zip
    my @files = ($array . ".na"    . $vers . ".annot.csv",
                 $array . ".cn.na" . $vers . ".annot.csv");
    my $subUrl = $baseUrl . "/genotyping/";
    my %urls   = map { $_ => $subUrl.$_.".zip" } @files;
    
    my $type  = "Map";
    my @ns    = ("EntrezGene", "AffyProbeSet");
    my %fbits = (type => $type,  mod  => $arrayToken,
                 ns1  => $ns[0],  ns2 => $ns[1],
                 auth => $auth,  vers => $versToken);

    my $byCtx = $args->{bycontext} ? 1 : 0;

    my @urls;
    my $meta = {
        MatrixType => $fbits{type},
        Modifier   => $fbits{mod},
        Version    => $fbits{vers},
        Source     => [values %urls],
    };
    my $fmeta = { %{$stashMeta}, %{$meta} };
    my $trg   = &primary_path(%fbits);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing $ns[0] $fbits{mod} map:", $trg);
        &post_process( %fbits, meta => $fmeta );
        return $trg;
    }

    my ($lvls, $scH) = &_factorize_levels('contextlevels');

    my ($rnum, $cnum, $nznum) = (0,0,0);
    my (%rmeta, %cmeta, %species);
    
    foreach my $file (@files) {
        my $csv  = &local_file( $file );
        my $zUrl = $urls{$file};

        my $zip  = &get_url($zUrl);
        if (&source_needs_recovery($csv)) {
            unzip $zip => $csv;
        }
        warn "  Parsing $csv\n";
        my $inData = 0;
        my (%tags, $header, $idCol, $gCol);
        my $io = Text::CSV->new({ sep_char => ',' });
        open(CSV, "<$csv") || die "Failed to read CSV\n  $csv  \n$!\n ";
        my $n = 0;
        while (<CSV>) {
            s/[\n\r]+$//;
            if (!$header) {
                if (/^</) {
                    &_unregistered_error( $zip, $csv, $zUrl );
                    last;
                } elsif (/^#%(.+?)=(.+)\s*$/) {
                    $tags{$1} = $2;
                } elsif (/^"/) {
                    my $stat = $io->parse($_);
                    if ($stat == 1) {
                        $header = [ $io->fields() ];
                        my %lu = map { lc($header->[$_]) => 
                                           $_ + 1 } (0..$#{$header});
                        $idCol = $lu{"probe set id"};
                        $gCol  = $lu{"associated gene"};
                        unless ($idCol && $gCol) {
                            warn join("\n    ", "Failed to identify probe and gene metadata columns. Observed columns:", @{$header}, "\n");
                            last;
                        }
                        $idCol--;
                        $gCol--;
                    } else {
                        die "Error parsing presumptive header:\n  $_\n  ";
                    }
                }
                next;
            }
            ## If here, we are processing the data block in the CSV
            my $stat = $io->parse($_);
            if ($stat != 1) {
                warn "Failed to parse row: $_\n";
                next;
            }
            my @row = $io->fields();
            my $id = $row[ $idCol ];
            next unless ($id);
            ## Columns will be probes / feature - there are about a
            ## million of them, so we're not going to attempt to bulk
            ## up the file with extra probe-level metadata.
            my $cm = $cmeta{$id} ||= {
                name  => $id,
                order => ++$cnum,
            };
            my $cn   = $cm->{order};
            foreach my $gbit (split(' /// ', $row[ $gCol ])) {
                ## Transcript ID, Context, Distance, UniGene, Symbol,
                ## EntrezGene ID, Description
                my ($rna, $ctx, $dist, $UG, $sym, $gid, $desc) =
                    split(' // ', $gbit);
                next unless ($gid && $gid =~ /^\d+$/);
                my $rm = $rmeta{$gid} ||= {
                    ## Build the row metadata hash
                    name        => $gid,
                    order       => ++$rnum,
                    Symbol      => $sym  eq '---' ? '' : $sym,
                    Description => $desc eq '---' ? '' : $desc,
                    ## ... plus row-to-col hash:
                    hits        => {},
                };
            
                my $hits = $rm->{hits};
                if ($byCtx) {
                    my $sc   = $scH->{uc($ctx)} || 1;
                    warn "[!!] The context '$ctx' has not been ordered in the -contextlevels parameter, so it will be classed as '$lvls->[0]'\n" if ($sc == 1 && !$tasks{"UnkCtx"}{$ctx}++);
                    ## Score assignments by the factorized context
                    if (!$hits->{$cn} || $hits->{$cn} < $sc) {
                        ## Record if this is the first time we've seen the
                        ## connection, or update if we found a 'better'
                        ## context to represent it (eg was 'downstream'
                        ## before, but we just found an 'exon'
                        ## connection).
                        $hits->{$cn} = $sc;
                    }
                } else {
                    ## Score assignements by the distance from marker to gene
                    if ($dist == 0) {
                        ## Sparse matrices consider a value of zero to
                        ## be "not there", so we will remap these
                        ## values to token negative numbers.
                        if ($ctx eq 'intron') {
                            $dist = -1;
                        } else {
                            ## All non-intronic zero-length contexts
                            ## will have the "best" score of -2. This
                            ## includes splice-site, exon, UTR-3,
                            ## UTR-5, cds, synon, missense and
                            ## nonsense
                            $dist = -2;
                        }
                    }
                    if (!$hits->{$cn} || $hits->{$cn} > $dist) {
                        ## Record the shortest distance
                        $hits->{$cn} = $dist;
                    }
                }
            }
            ++$n;
            warn sprintf("   %10d: %d genes, %d probes\n", $n, $rnum, $cnum)
                unless ($n % 5000);
            last if ($args->{limit} && $n >= $args->{limit});
        }
        close CSV;
        if (my $gs = $tags{'genome-species'}) {
            $species{$gs}++;
        }
    }
    
    my @rowIds = map { $_->{name} } sort 
    { $a->{order} <=> $b->{order} } values %rmeta;
    my @colIds = map { $_->{name} } sort 
    { $a->{order} <=> $b->{order} } values %cmeta;

    ## Tally up factor counts and total non-zero counts
    my %counts;
    foreach my $dat (values %rmeta) {
        my @u   = values %{$dat->{hits}};
        $nznum += $#u + 1;
        if ($byCtx) {
            map { $counts{ $lvls->[$_-1]}++ } @u;
        }
    }

    my $tmp = "$trg.tmp";
    open(MTX, ">$tmp") || &death("Failed to write $ns[1] map", $tmp, $!);
    
    print MTX &_initial_mtx_block
        ( "Mapping", $rnum, $cnum, $nznum, "$auth $array $ns[0]-to-$ns[1] Map",
          "Accession conversion from $ns[0] to $array $ns[1]",
          $byCtx ? "Factor levels representing probe set location relative to gene" : "Genomic distance from probe to gene; -2 = exon, -1 = intron",
          $ns[0], $ns[1]);
    
    print MTX &_dim_block({
        %{$fmeta},
        RowDim    => $ns[0],
        RowUrl    => $nsUrl->{$ns[0]},
        ColDim    => $ns[1],
        ColUrl    => $nsUrl->{$ns[1]}, 
        Authority => $authLong });
    
    print MTX &_citation_MTX();
    print MTX &_species_MTX( [ sort keys %species ] );
    if ($byCtx) {
        print MTX &_filter_block({
            AutoFilterComment => "NOTE: The standard automatic filters for this matrix will suppress upstream / downstream relationships. Use the \$reset() method to restore all connections.",
            TossLevel => "[,][unknown,upstream,downstream] ## By default, only features 'inside' genes are mapped to loci" });
        print MTX &_factor_map_block
            ( $lvls, $scH, "Probe Context",
              ["Location of probe relative to exon structure of gene"],
              \%counts);
    } else {
        print MTX &_filter_block({
            AutoFilterComment => "NOTE: The standard automatic filters for this matrix will suppress probesets that are more than 10kb away from a gene. Use the \$reset() method to restore all connections.",
            MaxScore => "10000 ## Exclude any probes that are more than 10kb away from the gene", });
    }
    
    print MTX &_rowcol_meta_comment_block(\%defColDef);

    ## Note row metadata
    print MTX &_generic_meta_block(\@rowIds, 'Row', \%rmeta, \@stndMeta);

    ## Column metadata
    print MTX "% $bar\n";
    print MTX &_generic_meta_block(\@colIds, 'Col', \%cmeta, []);

    print MTX &_triple_header_block( $rnum, $cnum, $nznum );
    for my $i (0..$#rowIds) {
        while (my ($j, $sc) = each %{$rmeta{$rowIds[$i]}{hits}}) {
            printf(MTX "%d %d %d\n", $i+1, $j, $sc);
        }
    }
    close MTX;

    rename($tmp, $trg);
    &msg("Generated $ns[0] to $ns[1] mapping", $trg);

    &post_process( %fbits, meta => $fmeta );
    return $trg;
}

sub _parse_csv_row {
    my $row = shift;
    while ($row =~ /^,?"(.+?)"/) {
    }
}

sub _unregistered_error {
    my ($zip, $file, $url) = @_;
    warn "
===================== Likely Authentication Error =====================
It looks like this file is an HTML document:

  $file

This generally means that the file was automatically downloaded, but
NetAffx did not recognize you as a registered user. You will need to:

1. Delete the cached CSV file:
   $file

2. Register with NetAffx (if not done already) then log in:
   $naRegister

3. Download this file:
     $url
   and save to this location (replacing the existing file)
     $zip

4. Re-run this script

=======================================================================
";
    
}

sub _citation_MTX {
    return &_default_parameter( "Citation", "Annotation information recovered from Affymetrix technical support page: http://www.affymetrix.com/support/technical/byproduct.affx")."%\n";
}
