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

my $tinyScore  = 0.001;
my @stndMeta   = qw(Symbol Description);
my %defColDef = ( Symbol => "Entrez Gene symbol as provided by Affymetrix",
                  Description => "Short descriptive text for a gene or probe" );

my $nsUrl = {
    EntrezGene     => 'https://www.ncbi.nlm.nih.gov/gene/%s', # Integer IDs
};


my $aReq  = $args->{array}   || "";
my $vers  = $args->{version} || "";
my $limit = $args->{limit}   || 0;

$clobber ||= 1 if ($limit);

if ($vers =~ /(\d+)/) { $vers = $1; }

my ($array, $arrayToken);
if ($aReq =~ /(\S+?)(\.cn)?.na(\d+)/) {
    ## Explicit file specification
    $array = $1;
    warn "[!] You specified a version as '$vers', but it will be replaced with '$3'\n" if ($vers && $vers ne "$3");
    $vers  = $3;
} elsif ($aReq =~ /(gw.?6|genomewide.*6)/i) {
    # http://www.affymetrix.com/Auth/analysis/downloads/na34/genotyping/GenomeWideSNP_6.cn.na34.annot.csv.zip
    $array      = "GenomeWideSNP_6";
    $arrayToken = "GW6";
    $vers ||= 34; # Most recent (Sep 2013) as of March 2018
} elsif ($aReq =~ /133/) {
    $array = "HG-U133";
    if ($aReq =~ /a/) {
        $array .= "A";
    } elsif ($aReq =~ /b/i) {
        $array .= "B";
    }
    $array = "HT_$array" if ($aReq =~ /ht/i);
    $array .= "_Plus"  if ($aReq =~ /plus/i);
    if ($aReq =~ /pm/i) {
        $array .= "_PM";
    } elsif ($aReq =~ /2/) {
        $array .= "_2";
    }

    # ivt/HG-U133A.na35.annot.csv.zip
    # ivt/HG-U133A_2.na35.annot.csv.zip
    # ivt/HG-U133B.na35.annot.csv.zip
    # ivt/HG-U133_Plus_2.na35.annot.csv.zip
    # ivt/HT_HG-U133A.na35.annot.csv.zip
    # ivt/HT_HG-U133B.na35.annot.csv.zip
    # ivt/HT_HG-U133_Plus_PM.na35.annot.csv.zip
      
} elsif ($aReq) {
    warn "Array request '$aReq' not recognized\n";
}

$arrayToken ||= $array;


my $naRegister = "http://www.affymetrix.com/site/terms.affx?dest=register&buttons=on";

if (!$array || !$vers || $args->{h} || $args->{help}) {
    warn "
Usage:

This program will generate gene sets for the SetFisher enrichment
package. It will recover information from the Affymetrix web site
based on an array design and a NetAffx version that you provide.

Required Argument:

    -array The Affymetrix array design desired. Arrays currently supported:

      GenomeWideSNP_6
      HT_HG-U133_Plus_PM
      HG-U133_Plus_2
      HG-U133A_2
      HG-U133A
      HG-U133B

   -version The NetAffx version to use. If not provided, this will be
            automatically assigned to the most recent version *prior*
            to NA36. It appears that there is a new system in place by
            Thermo Fisher Scientific (Affymertix's parent organization
            as of writing), and I don't have time at the moment to
            explore how it differs from the 'legacy' system used by
            NetAffx.

Optional Arguments:

   -tmpdir Default '$defTmp'. Directory holding downloaded files

      -dir Default '.'.  Output directory that will contain generated
           matrix and metadata files

  -clobber Default 0, which will preserve any already-generated
           files. Set to 1 to regenerate matrix files, and any value
           greater than 1 to also re-download base files from Entrez.

    -limit Default 0. If a non-zero value, will only process that
           number of probes. Will set -clobber to 1, and will not
           do post-processing (loading to stash, if available).

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

NOTE: Access to the annotation files requires registering with
NetAffx. You can review the Affymetrix terms of use here, as well as
register a new account.

  $naRegister

The script will attempt to download the appropriate files for you, but
will fail in ~100% of the cases. It will detect this failure, however,
and provide instructions on how to download the files manually (as
well as how you can register an account).

";

    warn "
You will need to specify an array. You can find a list of arrays here:

    http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays

You will need to match the format of the array name to that in the
link pointing to the annotation file. For example, if the link to the
array is 

   .../ivt/HG-U133A_2.na35.annot.csv.zip

... then the array would be 'HG-U133A_2'

    " unless ($array);

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


&msg("Array: $array");
&msg("NetAffx Version: $vers");
&msg("Output directory:",  $outDir);


my $baseUrl = "http://www.affymetrix.com/Auth/analysis/downloads/na" . $vers;


if ($arrayToken eq 'GW6') {
    $nsUrl->{AffyProbeSet} = "https://www.affymetrix.com/analysis/netaffx/mappingfullrecord.affx?pk=${array}:%s";
    &parse_genotyping();
} else {
    $nsUrl->{AffyProbeSet} = "https://www.affymetrix.com/analysis/netaffx/fullrecord.affx?pk=${array}:%s";
   &parse_ivt();
}

sub parse_ivt {
    my $subUrl = $baseUrl . "/ivt/";
    ## Should be just one annotation file for these arrays
    my $file = $array . ".na"    . $vers . ".annot.csv";
    my $url  = $subUrl . $file . ".zip";

    my $type  = "Map";
    my $ns2   = "AffyProbeSet";

    ## To calculate the fraction of probes matching a target RNA, we
    ## need the total number of probes in a probeset
    my $psCount = &_probe_set_counts( $array );


    ## We will simultaneously parse RefSeq and Ensembl transcripts
    ## into separate files
    my %mtxInfo;
    my $cnum = 0;
    my $sharedCmeta = {};
    my $nsLU = { refseq  => "RefSeqRNA",
                 ensembl => "EnsemblRNA" };
    foreach my $ns1 ("RefSeqRNA", "EntrezGene", "EnsemblRNA", "EnsemblGene") {
        my %fbits = (type => $type,  mod  => $arrayToken,
                     ns1  => $ns1,  ns2 => $ns2,
                     auth => $auth,  vers => $versToken);

        my $meta = {
            MatrixType => $fbits{type},
            Modifier   => $fbits{mod},
            Version    => $fbits{vers},
            Source     => $url,
        };
        my $fmeta = { %{$stashMeta}, %{$meta} };
        my $trg   = &primary_path(%fbits);
        if ($limit) {
            $trg =~ s/\.mtx$/-LIMIT$limit.mtx/;
        }
        unless (&output_needs_creation($trg)) {
            &msg("Keeping existing $ns1 $fbits{mod} map:", $trg);
            &post_process( %fbits, meta => $fmeta ) unless ($limit);
            next;
        }
        $mtxInfo{$ns1} = {
            fbits => \%fbits,
            fmeta => $fmeta,
            trg   => $trg,
            rnum  => 0,
            rmeta => {},
        };
    }

    ### MAPPING TO GENES ###

    ## TL;DR: No gene linkages available for EntrezGene

    ## If we wish to associate fractional probe alignment with each
    ## assignment (that is, #probesMatched / #probesInSet) we do not
    ## have reliable data to do so at the gene/locus level. Matched
    ## probe counts are provided for each transcript summarized in the
    ## "Transcript Assignments" column. It appears that for Ensembl
    ## the ENSG gene ID associated with the ENST transcript ID is
    ## available. However, the RefSeq IDs do not contain similar
    ## connections to the EntrezGene accessions.
    
    
    #my ($rnum, $cnum, $nznum) = (0,0,0);
    #    my (%rmeta, %cmeta);

    my $csv  = &local_file( $file );

    my $zip  = &get_url($url);
    if (&source_needs_recovery($csv)) {
        unzip $zip => $csv;
    }
    warn "  Parsing $csv\n";
    my $info = {
        zip    => $zip,
        csv    => $csv,
        url    => $url,
        tags   => {},
        header => "",
        ## Column indices:
        idCol  => undef, # Probe ID; should be unique
        syCol  => undef, # Probeset-level symbol; 1D list
        trgCol => undef, # Target description; should be unique?
        taCol  => undef, # Transcript annotation; 2D array
        egCol  => undef, # Entrez Gene ID; 1D list
        clCol  => undef, # Chromosome location; 2D array
        gtCol  => undef, # Gene Title, presumed 1D list
    };

    my @colMeta = ();
    my %colDefs = %defColDef;
    
    my ($idCol, $trgCol, $syCol, $egCol, $taCol, $clCol, $gtCol);
    $info->{checkhead} = sub {
        my $info   = shift;
        my $header = $info->{header};
        my %lu     = map { lc($header->[$_]) => 
                               $_ + 1 } (0..$#{$header});
        my $ic = $lu{"probe set id"};
        my $tc = $lu{"transcript assignments"};
        
        unless ($ic && $tc) {
            warn join("\n    ", "Failed to identify probe and transcript metadata columns. Observed columns:", @{$header}, "\n");
            return "Unknown Columns";
        }
        $idCol = $info->{idCol} = --$ic;
        $taCol = $info->{taCol} = --$tc;
        if (my $syc = $lu{"gene symbol"}) {
            $syCol = $info->{syCol} = --$syc;
        } else {
            warn "  [-] Note: Could not identify Symbol column\n";
        }
        if (my $gtc = $lu{"gene title"}) {
            $gtCol = $info->{gtCol} = --$gtc;
        }
        if (my $trc = $lu{"target description"}) {
            $trgCol = $info->{trgCol} = --$trc;
        } else {
            warn "  [-] Note: Could not identify Target Description column\n";
        }
        if (my $egc = $lu{"entrez gene"}) {
            $egCol = $info->{egCol} = --$egc;
            push @colMeta, "GeneCount";
            $colDefs{GeneCount} = "Number of distinct Entrez genes reported by Affymetrix. Will be equal to GenomeCount in nearly all cases.";
        } else {
            warn "  [-] Note: Could not identify Entrez Gene column\n";
        }
        if (my $clc = $lu{"chromosomal location"}) {
            $clCol = $info->{clCol} = --$clc;
            push @colMeta, "GenomeCount";
            $colDefs{GenomeCount} = "Number of distinct chromsomal locations reported by Affymetrix. Will be equal to GeneCount in nearly all cases.";
        } else {
            warn "  [-] Note: Could not identify Chromosomal Location column\n";
        }
        push @colMeta, ("Symbol", "Description");
        return 0;
    };


    my %noDescTags;
    
    my $io = Text::CSV->new({ sep_char => ',' });
    open(CSV, "<$csv") || die "Failed to read CSV\n  $csv  \n$!\n ";
    my $n = 0;
    while (<CSV>) {
        s/[\n\r]+$//;
        if (!$info->{header}) {
            if (my $issue = &_parse_csv_header( $_, $info, $io )) {
                warn "[!!] File not parsed: $issue\n";
                last;
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

        ## All four matrices share the same probe set information:
        my $desc = "";
        if (defined $trgCol) {
            $desc = $row[ $trgCol ] || "";
            ## See if we can tidy up the description a bit
            my $ft = &_parse_fasta_tags( $desc );
            if (my $def = $ft->{def} || $ft->{definition}) {
                ## If /DEF or /DEFINITION are defined, use that, it
                ## generally seems to be the most reasonable
                ## description.
                $desc = $def;
            } else {
                my $res = $ft->{RESIDUAL}; # text after tags removed
                my $ug  = $ft->{ug} || ""; # Unigene text
                my $ugXtra = $ug; $ugXtra =~ s/^[a-z]+.\d+\s*//i;
                if (length($ugXtra) > 12) {
                    ## Looks like we have at least a modest UniGene
                    ## description (more than just accession). Use it:
                    $desc = $ug;
                } elsif ($res) {
                    ## Otherwise, use the residual text
                    $desc = $res;
               }
            }
        }
        $desc =~ s/\s*\.\s*$//; # Trailing whitespace / period
        ## In some cases, the description is not terribly helpful -
        ## just an accession number, or a sterile description of an
        ## EST cluster. In these cases, let's take the "Gene Title"
        ## information if it's available
        if (($desc =~ /^[a-z0-9_\.]+$/i || $desc =~ /^Cluster Incl. /)
            && defined $gtCol) {
            my @gt = &_2d_array($row[$gtCol]);
            if ($gt[0][0]) {
                my $better = $gt[0][0];
                ## More than one? just note how many more to avoid
                ## spammy descriptions eg (11715189_s_at):
                ## H3 histone, family 3A / H3 histone, family 3B (H3.3B) / H3 histone, family 3C / microRNA 4738
                my $xtra = $#gt;
                $better .= sprintf(" - plus %d other gene%s", $xtra, $xtra == 1 ? '': 's') if ($xtra);
                $desc = "$better ($desc)";
            }
        }
 
        my $chrCnt = "";
        if (defined $clCol) {
            my @cl =  &_2d_array($row[ $syCol ]);
            $chrCnt = $#cl + 1;
        }

        my $cm = $sharedCmeta->{$id} ||= {
            name   => $id,
            order  => ++$cnum,
            Symbol => defined $syCol ? &_2d_array($row[ $syCol ], ',') : "",
            GenomeCount => $chrCnt,
            Description => $desc,
        };
        my $cn     = $cm->{order};
        ## Total number of probes for this probeset:
        my $prbCnt = $psCount->{$id} || 0;
        
        if (defined $egCol && $mtxInfo{EntrezGene}) {
            ## Associate with EntrezGene
            my $rmeta = $mtxInfo{EntrezGene}{rmeta};
            $cm->{GeneCount} = 0;
            foreach my $egd (&_2d_array($row[ $egCol ])) {
                my $gid = $egd->[0];
                if ($gid && $gid =~ /^\d+$/) {
                    ## Keep a tally of the distinct genes mapped to probeset
                    $cm->{GeneCount}++;
                    my $rm = $rmeta->{$gid} ||= {
                        ## Build the row metadata hash
                        name        => $gid,
                        order       => ++$mtxInfo{EntrezGene}{rnum},
                        ## ... plus row-to-col hash:
                        hits        => {},
                    };
                    ## We can't associate EntrezGene with a score, so
                    ## just assign all connections a value of -1; A
                    ## negative number is used to note "Connection
                    ## reported, unknown similarity"
                    $rm->{hits}{$cn} = -1;
                }
            }
        }
        
        ## Check each transcript associated with the probeset:
        foreach my $taDat (&_2d_array($row[ $taCol ])) {
            if ($taDat->[2] eq 'refseq' && $mtxInfo{RefSeqRNA}) {
                ## Associate probeset with RefSeq RNA
                my $rid = $taDat->[0];
                if ($rid =~ /^[NX][MR]_(\d{6}|\d{9})$/) {
                    ## Looks ok
                    unless ($mtxInfo{RefSeqRNA}{rmeta}{$rid}) {
                        ## Try to sneak some structured info out of
                        ## the description.
                        my $desc = $taDat->[1] || "";
                        $desc =~ s/\s*\.\s*$//; # don't need '.' at end
                        my $sym = "";
                        ## We're trying to pull out the
                        ## frequently-embedded RefSeq gene symbol, eg
                        ## 'YY1' in this description:
                        ##     YY1 transcription factor (YY1), mRNA.
                        
                        ## This is messy because many RNAs have
                        ## additional parenthetical short strings
                        ## followed by a comma and space. So we will
                        ## look for specific locus flags after the
                        ## comma to identify the symbol notation:
                        if ($desc =~ / \(([^\)]+)\), (transcript|mRNA|non-coding|long non-coding|misc_RNA|ncRNA|microRNA|antisense RNA|partial mRNA|small nucleolar RNA|guide RNA|ribosomal RNA|RNase P RNA|partial misc_RNA)/) {
                            $sym = $1;
                        }
                        ## Examples of problem notation:
                        ##   ... sub-family F (GCN20), member 3 (ABCF3), mRNA
                        ##   ... A (SII), 2 (TCEA2), transcript variant 1, mRNA
                        my $rm = $mtxInfo{RefSeqRNA}{rmeta}{$rid} ||= {
                            ## Build the row metadata hash
                            name        => $rid,
                            order       => ++$mtxInfo{RefSeqRNA}{rnum},
                            Symbol      => $sym,
                            Description => $desc,
                            ## ... plus row-to-col hash:
                            hits        => {},
                        };
                    }
                    my $rm = $mtxInfo{RefSeqRNA}{rmeta}{$rid};
                    my $hits = $rm->{hits};
                    my $sc   = &_probe_score($taDat->[3], $prbCnt);
                    if (!$hits->{$cn} || $hits->{$cn} < $sc) {
                        ## Record the best score observed for this transcript
                        $hits->{$cn} = $sc;
                    }
                }
            } elsif ($taDat->[2] eq 'ensembl') {
                ## Ensembl transcript
                my $rid = $taDat->[0];
                if ($rid =~ /^ENS[A-Z]{0,3}T\d+$/) {
                    ## Looks like an ok Ensembl RNA ID
                    my $sc   = &_probe_score($taDat->[3], $prbCnt);
                    if ($mtxInfo{EnsemblRNA}) {
                        ## We are building Ensembl RNA matrix
                        my $rm = $mtxInfo{EnsemblRNA}{rmeta}{$rid} ||= {
                            ## Build the row metadata hash
                            name        => $rid,
                            order       => ++$mtxInfo{EnsemblRNA}{rnum},
                            ## The description field is NOT useful for ENST!
                            ## ... plus row-to-col hash:
                            hits        => {},
                        };
                        my $hits = $rm->{hits};
                        if (!$hits->{$cn} || $hits->{$cn} < $sc) {
                            ## Record the best score observed
                            $hits->{$cn} = $sc;
                        }
                    }
                    if ($mtxInfo{EnsemblGene}) {
                        ## We are building Ensembl Gene matrix; The
                        ## description was not useful for annotating
                        ## the transcripts, but it does seem to
                        ## contain gene accessions. See if we can find one
                        if ($taDat->[1] =~ /gene:(ENS[A-Z]{0,3}G\d+)/) {
                            my $gid = $1;
                            my $rm = $mtxInfo{EnsemblGene}{rmeta}{$gid} ||= {
                                ## Build the row metadata hash
                                name        => $gid,
                                order       => ++$mtxInfo{EnsemblGene}{rnum},
                                ## The description field is NOT useful for ENST!
                                ## ... plus row-to-col hash:
                                hits        => {},
                            };
                            my $hits = $rm->{hits};
                            if (!$hits->{$cn} || $hits->{$cn} < $sc) {
                                ## Record the best score observed;
                                ## This is probably the only class of
                                ## matrix where explicitly selecting
                                ## the best score is relevant, since
                                ## multiple transcripts might result
                                ## in a variety of scores associated
                                ## with the gene.
                                $hits->{$cn} = $sc;
                            }
                        }
                    }
                }
            }
        }
        ++$n;
        warn sprintf("   %10d: %s\n", $n, $id) unless ($n % 1000);
        last if ($limit && $n >= $limit);
    }
    close CSV;

    my @ndt = sort { $noDescTags{$b} <=> $noDescTags{$a} } keys %noDescTags;
    if ($#ndt != -1) {
        warn join("\n", "Fasta tags on targets that didn't have obvious description tags:", map { sprintf("    %5d %s", $noDescTags{$_}, $_) } @ndt)."\n";
    }
    

    foreach my $ns1 (sort keys %mtxInfo) {
        my $mti  = $mtxInfo{$ns1};
        my $rnum = $mti->{rnum};
        if ($rnum == 0) {
            warn "  [??] No connections found to $ns1\n";
            next;
        }
        my $scDesc = $ns1 eq 'EntreGene' ? "All scores are -1" : "Fraction of probes in probeset matching target; -1 indicates unknown, $tinyScore indicates explicitly reported zero probes matching";

        ## Only RefSeq has reliable row metadata:
        my $rmCols = $ns1 eq 'RefSeqRNA' ? ['Symbol', 'Description'] : [];

        my $filtCB = $ns1 eq 'EntrezGene' ?
            sub {
                ## Nothing to do with EntrezGene, we don't have
                ## quantitative information to make decent filters.
        } : sub {
            my ($fh, $counts) = @_;
            print $fh &_filter_block({
                AutoFilterComment => "NOTE: The standard automatic filters for this matrix will only preserve connections where at least 80% of the probes matching the target. To recover probests with weaker connections (or unquantified ones), use the \$reset() method to restore all connections.",
                MinScore => "0.8 ## Only include probesets that have at least 80% of probes matching the target RNA", });            
        };

        &_process_data($mti->{rmeta}, $sharedCmeta, $rnum, $cnum, 0,
                       $rmCols, \@colMeta,
                       $mti->{trg}, $mti->{fbits}, $mti->{fmeta}, 
                       $scDesc, $filtCB,
                       $info, \%colDefs);
    }
}

sub parse_genotyping {
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
    $trg =~ s/\.mtx$/-LIMIT$limit.mtx/ if ($limit);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing $ns[0] $fbits{mod} map:", $trg);
        &post_process( %fbits, meta => $fmeta ) unless ($limit);
        return $trg;
    }

    my ($lvls, $scH) = &_factorize_levels('contextlevels');

    my ($rnum, $cnum) = (0,0);
    my (%rmeta, %cmeta, %topInfo);
    
    foreach my $file (@files) {
        my $csv  = &local_file( $file );
        my $zUrl = $urls{$file};

        my $zip  = &get_url($zUrl);
        if (&source_needs_recovery($csv)) {
            unzip $zip => $csv;
        }
        warn "  Parsing $csv\n";
        my $info = {
            zip    => $zip,
            csv    => $csv,
            url    => $zUrl,
            tags   => {},
            header => "",
            idCol  => 0,
            gCol   => 0,
        };
        $info->{checkhead} = sub {
            my $info   = shift;
            my $header = $info->{header};
            my %lu     = map { lc($header->[$_]) => 
                                   $_ + 1 } (0..$#{$header});
            my $idCol = $lu{"probe set id"};
            my $gCol  = $lu{"associated gene"};
            unless ($idCol && $gCol) {
                warn join("\n    ", "Failed to identify probe and gene metadata columns. Observed columns:", @{$header}, "\n");
                return "Unknown Columns";
            }
            $info->{idCol} = $idCol--;
            $info->{gCol}  = $gCol--;
            ## Aggregate the tags from the two files
            my %aggInfo;
            while (my ($k, $v) = each %{$info->{tags}}) {
                $aggInfo{$k}{$v} = 1 if ($k && $v);
            }
            while (my ($k, $vH) = each %aggInfo) {
                my @u = keys %{$vH};
                $topInfo{$k} = $#u == 0 ? $u[1] : \@u;
            }
            return 0;
        };
        my ($idCol, $gCol);
        my $io = Text::CSV->new({ sep_char => ',' });
        open(CSV, "<$csv") || die "Failed to read CSV\n  $csv  \n$!\n ";
        my $n = 0;
        while (<CSV>) {
            s/[\n\r]+$//;
            if (!$info->{header}) {
                if (my $issue = &_parse_csv_header( $_, $info, $io )) {
                    warn "[!!] File not parsed: $issue\n";
                    last;
                }
                $idCol ||= $info->{idCol};
                $gCol  ||= $info->{gCol};
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
            ## million of them in each file (935k SNPs, 946k CNVs), so
            ## we're not going to attempt to bulk up the file with
            ## extra probe-level metadata.
            my $cm = $cmeta{$id} ||= {
                name  => $id,
                order => ++$cnum,
            };
            my $cn   = $cm->{order};
            foreach my $gdat (&_2d_array($row[ $gCol ])) {
                ## Transcript ID, Context, Distance, UniGene, Symbol,
                ## EntrezGene ID, Description
                my ($rna, $ctx, $dist, $UG, $sym, $gid, $desc) = @{$gdat};
                next unless ($gid && $gid =~ /^\d+$/);
                my $rm = $rmeta{$gid} ||= {
                    ## Build the row metadata hash
                    name        => $gid,
                    order       => ++$rnum,
                    Symbol      => $sym,
                    Description => $desc,
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
            last if ($limit && $n >= $limit);
        }
        close CSV;
    }
    my $scDesc = $byCtx ? "Factor levels representing probe set location relative to gene" : "Genomic distance from probe to gene; -2 = exon, -1 = intron";
    my $filtCB = $byCtx ?
        sub {
            my ($fh, $counts) = @_;
            print $fh &_filter_block({
                AutoFilterComment => "NOTE: The standard automatic filters for this matrix will suppress upstream / downstream relationships. Use the \$reset() method to restore all connections.",
                TossLevel => "[,][unknown,upstream,downstream] ## By default, only features 'inside' genes are mapped to loci" });
            print $fh &_factor_map_block
                ( $lvls, $scH, "Probe Context",
                  ["Location of probe relative to exon structure of gene"],
                  $counts);
    } : sub {
        my ($fh, $counts) = @_;
        print $fh &_filter_block({
            AutoFilterComment => "NOTE: The standard automatic filters for this matrix will suppress probesets that are more than 10kb away from a gene. Use the \$reset() method to restore all connections.",
            MaxScore => "10000 ## Exclude any probes that are more than 10kb away from the gene", });
    };
    &_process_data(\%rmeta, \%cmeta, $rnum, $cnum, $byCtx ? 'contextlevels': 0,
                   \@stndMeta, [], $trg, \%fbits, $fmeta, $scDesc, $filtCB,
                   \%topInfo);
}

sub _process_data {
    ## Row data, Col data, num rows, num cols, factor argument name,
    ## Row metadata columns, Col metadata columns, Output file target, 
    ## FileName bits, File metdata, Score descriptive text
    ## Filter block callback function, Parsing information
    my ($rmeta, $cmeta, $rnum, $cnum, $byFactor,
        $rmCols, $cmCols, $trg,
        $fbits, $fmeta, $scDesc,
        $filtCB, $info, $colDef) = @_;

    $colDef ||= \%defColDef;

    my @rowIds = map { $_->{name} } sort 
    { $a->{order} <=> $b->{order} } values %{$rmeta};
    my @colIds = map { $_->{name} } sort 
    { $a->{order} <=> $b->{order} } values %{$cmeta};

    my ($lvls, $scH);
    ($lvls, $scH) = &_factorize_levels($byFactor) if ($byFactor);

    ## Tally up factor counts and total non-zero counts
    my $nznum = 0;
    my %counts;
    foreach my $dat (values %{$rmeta}) {
        my @u   = values %{$dat->{hits}};
        $nznum += $#u + 1;
        if ($byFactor) {
            map { $counts{ $lvls->[$_-1]}++ } @u;
        }
    }

    my ($ns1, $ns2) = map { $fbits->{$_} } qw(ns1 ns2);

    my $tmp = "$trg.tmp";
    my $mtx;
    open($mtx, ">$tmp") || &death("Failed to write $ns2 map", $tmp, $!);
    
    print $mtx &_initial_mtx_block
        ( "Mapping", $rnum, $cnum, $nznum, "$auth $array $ns1-to-$ns2 Map",
          "Accession conversion from $ns1 to $array $ns2", $scDesc,
          $ns1, $ns2);
    
    print $mtx &_dim_block({
        %{$fmeta},
        RowDim    => $ns1,
        RowUrl    => $nsUrl->{$ns1},
        ColDim    => $ns2,
        ColUrl    => $nsUrl->{$ns2}, 
        Authority => $authLong });

    print $mtx &_citation_MTX();
    print $mtx &_species_MTX( $info->{tags}{'genome-species'} );
    &{$filtCB}( $mtx, \%counts );
    
    print $mtx &_rowcol_meta_comment_block($colDef);

    ## Note row metadata
    print $mtx &_generic_meta_block(\@rowIds, 'Row', $rmeta, $rmCols);

    ## Column metadata
    print $mtx "% $bar\n";
    print $mtx &_generic_meta_block(\@colIds, 'Col', $cmeta, $cmCols);

    print $mtx &_triple_header_block( $rnum, $cnum, $nznum );
    for my $i (0..$#rowIds) {
        while (my ($j, $sc) = each %{$rmeta->{$rowIds[$i]}{hits}}) {
            printf($mtx "%d %d %s\n", $i+1, $j, $sc);
        }
    }
    close $mtx;

    rename($tmp, $trg);
    &msg("Generated $ns1 to $ns2 mapping", $trg);

    &post_process( %{$fbits}, meta => $fmeta ) unless ($limit);
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

sub _parse_csv_header {
    my ($line, $info, $io) = @_;
    if (/^</) {
        
        &_unregistered_error( $info->{zip}, $info->{csv}, $info->{url} );
        return "Unregistered"
    } elsif (/^#%(.+?)=(.+)\s*$/) {
        ## Key-value tags at top of file
        $info->{tags}{$1} = $2;
    } elsif (/^"/) {
        my $stat = $io->parse($_);
        if ($stat == 1) {
            $info->{header} = [ $io->fields() ];
            return &{$info->{checkhead}}( $info );
            
        } else {
            die "Error parsing presumptive header:\n  $_\n  ";
        }
    }
    return 0;
}

sub _2d_array {
    ## Splits content of the form:
    ##    IPR001523 // Paired domain // 1.4E-75 /// IPR001523 // Paired domain // 1.1E-75 /// etc
    my ($val, $flatten, $subFlatten) = @_;
    if (!$val || $val eq '---') {
        return $flatten ? "" : ();
    }
    my @rv;
    foreach my $bit (split(' /// ', $val)) {
        push @rv, [ map { $_ eq '---' ? '' : $_ } split(' // ', $bit) ];
    }
    if ($flatten) {
        my $rv = join($flatten, map { join($subFlatten || ',', @{$_}) } @rv);
        return defined $rv ? $rv : "";
    } else {
        return @rv;
    }
}

sub _probe_set_counts {
    my $arr = shift;
    ## http://www.affymetrix.com/Auth/analysis/downloads/data/HG-U133A.probe_tab.zip
    my $file = $arr . ".probe_tab";
    my $url  = "http://www.affymetrix.com/Auth/analysis/downloads/data/".
        $file.".zip";
    my $tsv  = &local_file( $file );
    
    my $zip  = &get_url($url);
    if (&source_needs_recovery($tsv)) {
        unzip $zip => $tsv;
    }
    warn "  Counting probes in $arr probesets\n";
    open(PSEQ, "<$tsv") || die "Failed to read probe file\n  $tsv  \n$!\n ";
    my $head = <PSEQ>;
    if ($head =~ /^</) {
        &_unregistered_error( $zip, $tsv, $url );
        close PSEQ;
        die;
    }
    my %counts;
    while (<PSEQ>) {
        ## We just want to know how many entries (probes) are present
        ## for each probe set, so simply tally occurances of first
        ## column.
        my @bits = split(/\t/);
        $counts{$bits[0]}++;
    }
    close PSEQ;
    return \%counts;
}

sub _probe_score {
    ## Return a fractional score relating the number of probes
    ## presumed to hit a transcript to the fraction that actually
    ## do. So if a probeset has 13 probes, an 11 are found to actually
    ## align to the transcript, the score will be 11/13 = 0.846
    my ($probeCount, $totalCount) = @_;
    if (!$totalCount) {
        ## If we don't know how many probes the probeset has, use -1
        ## to represent "noted hit but unknown similarity"
        return -1;
    } elsif (!$probeCount) {
        ## The sparse matrix will exclude values of zero. So record
        ## zeros as a "tiny value" to allow explicitly zero hits to
        ## still be noted:
        return $tinyScore;
        ## ... it's often handy to know that a transcript was
        ## inspected but found to NOT match the probeset.
    } else {
        ## Otherwise, return the fraction of total probes matched, 3 sigfig
        return int(0.5 + 1000 * $probeCount / $totalCount) / 1000;
    }
}

sub _parse_fasta_tags {
    my $val = $_[0] || "";
    my %rv;
    ## These tags are kinda tough to identify and cleanly pull out,
    ## particularly when some entries are a combination of /TAG= and
    ## un-tagged free-form text. We'll find the end of the tag by
    ## looking for the start of the next tag, or the end of the
    ## string.
    while ($val =~ /(\/([a-z_]+)=(.*?))( \/([a-z_]+)=|\s*$)/i) {
        my ($rep, $k, $v) = ($1, $2, $3);
        $val =~ s/\Q$rep\E//; ## Excise the block we just found
        $rv{lc($k)} = $v;
    }
    $val =~ s/\s+$//;
    $rv{RESIDUAL} = $val if ($val); ## What is left over
    # die "$_[0] = ".join(" + ", keys %rv);
    return \%rv;
}
