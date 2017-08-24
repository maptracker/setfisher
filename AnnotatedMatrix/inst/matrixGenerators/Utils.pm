use strict;
use Data::Dumper;
use File::Basename;
use File::Path 'mkpath';
use POSIX 'strftime';
use File::Listing 'parse_dir';
use IO::Uncompress::Gunzip;
use Net::FTP;
use LWP::UserAgent;

our ($defaultArgs, $defTmp, $ua);
my $codeDir     = dirname($0);
our $bar        = "- " x 20;
our $mtxSep     = " :: ";


our $args = &parseargs({
    ## Tool-specific defaults:
    %{$defaultArgs},
    ## Defaults common to tools:
    tmpdir   => $defTmp,
    clobber  => 0,
    verify_hostname => 0,  });

my  $stash    = $args->{stash};
my  %tasks    = ();
our $tmpDir   = $args->{tmpdir}; $tmpDir =~ s/\/+$//;
our $clobber  = $args->{clobber} || 0;
our $maxAbst  = 150;
our $outDir   = $args->{dir} || "";
$outDir =~ s/\/+$//;

&mkpath([$tmpDir]) if ($tmpDir);

# die Dumper($args);

sub parseargs {
    ## Argument parsing
    my %rv;
    ## Normalize the defaults
    while (my ($k,$v) = each %{ $_[0] || {} }) {
        $k =~ s/^\-//;
        $rv{lc($k)} = $v if ($k && defined $v);
    }

    # Configuration files, both in the generators folder and the
    # user's home folder
    my @cFiles = map { sprintf("%s/.matrixGenerators.conf", $_) }
    ( $codeDir,  $ENV{HOME} || "UnSeThOmE");
    foreach my $cFile (@cFiles) {
        next unless (-s $cFile);
        ## A configuration file is present
        if (open(CF, "<$cFile")) {
            my @captured;
            while (<CF>) {
                next if (/^#/);
                s/[\n\r]+$//;
                if (/^\s*(.+?)\s*=\s*(.+?)\s*$/) {
                    my ($k, $v) = ($1, $2);
                    $k =~ s/^\-//;
                    if ($k && defined $v) {
                        $rv{lc($k)} = $v;
                        push @captured, $k;
                    }
                }
            }
            close CF;
            my $num = $#captured +1;
            &msg(sprintf("Loaded %d settings %sfrom %s", $num, 
                 $num ? '('.join(', ', @captured).') ' : '', $cFile));
        } else {
            &err("Failed to read configuration file '$cFile'", $!);
        }
    }
    
    ## Command line options
    my $i  = 0;
    while ($i <= $#ARGV) {
        my $key = lc($ARGV[$i]);
        $key =~ s/^\-+//;
        my $val = $i < $#ARGV && $ARGV[$i+1] !~ /^\-/ ? $ARGV[++$i] : 1;
        $rv{$key} = $val;
        $i++;
    }
    return \%rv;
}

sub msg {
    warn "[*] ".join("\n    ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}

sub err {
    warn "[!!] ERROR: ".join
        ("\n     ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}

sub death { &err(@_); die " -- "; }

sub _ftp {
    ## http://perlmeme.org/faqs/www/ftp_file_list.html
    my $site = shift || $args->{ftp};
    my $ftp = Net::FTP->new($site) ||
        &death("Failed to connect to FTP", $site, $!);
    my $pass = $site =~ /ncbi/ ? $args->{email} : "";
    $ftp->login("anonymous", $pass) ||
        &death("Failed to login to FTP", $site, $!);
    $ftp->binary();
    return $ftp;
}

=head2 SSL Issues with NCBI

At time of writing (9 May 2017) there was an issue with a mismatch
between the FTP site and the certificate:

   ERROR: certificate common name “*.ncbi.nlm.nih.gov” doesn’t match requested host name “ftp.ncbi.nih.gov”.

This will cause LWP to fail unless verify_hostname is set to false,
which is the current default. When NCBI supplies a certificate with
matching common name this option should be removed.

=cut

sub _ua {
    return $ua if ($ua);
    $ua       = LWP::UserAgent->new;
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
}

sub local_file {
    # Local temp file for a remote URL
    my ($uReq, $dReq) = @_;
    $dReq ||= basename($uReq);
    return "$tmpDir/$dReq";
}

sub ftp_url {
    my ($uReq, $domain) = @_;
    $domain ||= $args->{ftp};
    return join('/', "ftp:/", $domain, $uReq);
}

sub fetch_url {
    my ($uReq, $dReq) = @_;
    my $dest = &local_file( @_ );
    if (&source_needs_recovery($dest)) {
        ## File not yet downloaded, or request to re-download
        my $pending = 1;
        ## For NCBI it looks like the timeout is 60 seconds, and
        ## parsing is ~90 seconds per file on my system. Close and
        ## reinitialize to be assured of a 'live' connection:
        my $site;
        if ($uReq =~ /\/\/([^\/]+)\/(.+)/) {
            $site = $1;
            $uReq = $2;
        }
        my $ftp = &_ftp($site);
        my $tmp = "$dest.tmp";
        while ($pending) {
            $ftp->get($uReq, $tmp) || 
                &death("Failed to FTP file", $uReq, $ftp->message);
            if (-s $tmp) {
                rename($tmp, $dest);
                &msg("Downloaded $uReq:", $dest);
                if ($uReq =~ /(.+)\/([^\/]+)$/) {
                    ## Try to get remote timestamp
                    ## http://www.perlmonks.org/?node_id=493844
                    my ($dir, $file) = ($1, $2);
                    $ftp->cwd($dir);
                    my $ls = $ftp->dir();
                    foreach my $entry (parse_dir($ls)) {
                        my ($name, $type, $size, $mtime, $mode) = @$entry;
                        if ($name eq $file) {
                            my $time_string = strftime("%Y-%m-%d %H:%M:%S",
                                                       gmtime($mtime));
                            ## Change the date to match the FTP server
                            ## https://askubuntu.com/a/62496
                            system("touch --date='$time_string' '$dest'");
                            ## Eh. This is not being terribly precise,
                            ## possibly due to time zone issues,
                            ## possibly because HH::MM are not being
                            ## accurately returned. Should be within
                            ## 24hr of reality, which will work for
                            ## our needs.
                        }
                    }
                }
                undef $ftp;
                $pending = 0;
            } else {
                &err("Attempt $pending: Failed to recover remote file",
                     "Source: $uReq", "Destination: $dest", $ftp->message);
                return "" if ($pending++ >= 3);
                sleep(5);
            }
        }
    }
    return $dest;
}

sub get_url {
    my ($url, $dReq, $params) = @_;
    my $dest = &local_file( @_ );
    if (&source_needs_recovery($dest)) {
        &_ua();
        my $resp;
        if ($params) {
            $resp = $ua->post($url,':content_file' => $dest,Content => $params);
        } else {
            $resp = $ua->get($url, ':content_file' => $dest );
        }
        if ($resp->is_success) {
            if (-s $dest) {
                &msg("Downloaded $url :", $dest);
            } else {
                &death("Failed to recover URL:",
                       "   URL: $url",
                       "  File: $dest",
                       "Web request succesful, but file is empty");
            }
        } else {
            &death("Failed to download URL",
                   "   URL: $url",
                   "  File: $dest",
                   "Status: ".$resp->status_line);
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

my @mnThree    = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @mnFull     = qw(January February March April May June 
                    July August September October November December);
# These are not quite right, but close enough:
my $mnHash     = {
    winter       => 1,
    spring       => 4,
    spr          => 4,
    summer       => 7,
    sum          => 7,
    autumn       => 10,
    fall         => 10,
    easter       => 4,  # http://www.ncbi.nlm.nih.gov/pubmed/2019856
    christmas    => 12, # http://www.ncbi.nlm.nih.gov/pubmed/13954194
    '1st Quart'  => 1,  # http://www.ncbi.nlm.nih.gov/pubmed/10237212
    '2d Quart'   => 4,  # http://www.ncbi.nlm.nih.gov/pubmed/10248150
    '3d Quart'   => 7,  # http://www.ncbi.nlm.nih.gov/pubmed/10236507
    '4th Quart'  => 10, # http://www.ncbi.nlm.nih.gov/pubmed/10249456
    '-00 Winter' => 1,  # http://www.ncbi.nlm.nih.gov/pubmed/10711319
    '-94 Winter' => 1,  # http://www.ncbi.nlm.nih.gov/pubmed/11362190
    'N0v'        => 11, # http://www.ncbi.nlm.nih.gov/pubmed/5275180
};
foreach my $arr (\@mnThree, \@mnFull) {
    map { $mnHash->{lc($arr->[$_])} = $_ + 1 } (0..$#{$arr});
}
for my $i (1..12) {
    $mnHash->{$i} = $i;
    $mnHash->{sprintf("%02d", $i)} = $i;
}

sub parse_pubmed_date {
    my ($pmid, $y, $m, $d) = @_;
    my $date = "";
    if ($y && $#{$y} != -1) {
        &msg("Multiple years for PMID:$pmid") if ($#{$y} > 0);
        $y = $y->[0];
        if ($y =~ /^\d{4}$/) {
            $date = $y;
            if ($m && $#{$m} != -1) {
                &msg("Multiple months for PMID:$pmid") if ($#{$m} > 0);
                $m = $m->[0];
                if (my $std = $mnHash->{lc($m)}) {
                    $date = sprintf("%s-%02d", $date, $std);
                    if ($d && $#{$d} != -1) {
                        &msg("Multiple days for PMID:$pmid") if ($#{$d} > 0);
                        $d = $d->[0];
                        if ($d =~ /^\d{1,2}$/ && $d > 0 && $d < 32) {
                            $date = sprintf("%s-%02d", $date, $d);
                        } else {
                            &msg("Weird day for PMID:$pmid : '$d'");
                        }
                    }
                } else {
                    &msg("Weird month for PMID:$pmid : '$m'");
                }
            }
        } else {
            &msg("Weird year for PMID:$pmid : '$y'");
        }
    }
    return $date;
}

sub _triple_header_block {
    ## A comment header used just before triples are recorded in an
    ## MTX file, plus the row / col / non-zero count row
    my ($rnum, $cnum, $nznum) = @_;
return "% $bar
% Matrix triples : Row Col Score
% $bar
". sprintf("  %d %d %d\n", $rnum, $cnum, $nznum);
}

sub _mtx_comment_block {
    my $rv = "";
    my $wid = 78;
    for my $i (0..$#_) {
        my $com = $_[$i];
        if ($com eq "") {
            $rv .= "%\n";
        } elsif ($com =~ /^>(.+)<$/) {
            my $bit = $1;
            my $pad = int(($wid - length($bit))/2);
            $rv .= sprintf("%% %s%s\n", " " x $pad, $bit);
        } else {
            while (length($com) > $wid) {
                my $bit = substr($com, 0, $wid);
                $bit =~ s/\s+[^\s]+$//;
                $rv .= "% $bit\n";
                $com =~ s/^\Q$bit\E\s*//;
            }
            $rv .= "% $com\n" unless ($com eq "");
        }
    }
    return $rv;
}

sub _initial_mtx_block {
    my ($what, $rnum, $cnum, $nznum, $name, $desc, $scDesc,
        $rnm, $cnm, $rowAdj) = @_;
    ## Additional adjective describing row namespace (eg a species name)
    $rowAdj = $rowAdj ? "$rowAdj " : "";
    my $rv = "%%MatrixMarket matrix coordinate real general
% $what relating ${rowAdj}${rnm}s to ${cnm}s
% $rnum x $cnum sparse matrix with $nznum non-zero cells
% -- The 'AnnotatedMatrix' R package can parse these comments to decorate
%    the matrix:   myMatrix <- AnnotatedMatrix('path/to/this/file')
% Separator '$mtxSep'
%
";
    $rv .= &_default_parameter( "Name", $name );
    $rv .= &_default_parameter( "Description", $desc );
    $rv .= &_default_parameter( "ScoreDesc", $scDesc );
    $rv .="%\n";
    return $rv;
}

sub _dim_block {
    my $tags = shift || {};
    my $rv  = "";;
    my $rnm = $tags->{RowDim} || "";
    my $cnm = $tags->{ColDim} || "";
    $rv .= &_default_parameter( "RowDim", $rnm, "Rows are ${rnm}s" );
    $rv .= &_default_parameter( "RowUrl", $tags->{RowUrl});
    $rv .= &_default_parameter( "ColDim", $cnm, "Columns are ${cnm}s" );
    $rv .= &_default_parameter( "ColUrl", $tags->{ColUrl});
    $rv .= &_default_parameter( "Source", $tags->{Source},
                                "Data source(s):", " // ");
    $rv .= &_default_parameter( "Authority", $tags->{Authority});
    $rv .= "%\n" if ($rv);
    return $rv;
}


sub _rowcol_meta_comment_block {
    return "%% $bar
% Comment blocks for [Row Name] and [Col Name] follow (defining row
% and column names, plus metadata), followed finally by the triples
% that store the actual mappings.
%% $bar
";
}

sub gene_symbol_stats_block {
    my ($symbols, $ids, $lvls, $scH) = @_;

    my (@counts, %statuses, %oddChar);
    foreach my $sdat (values %{$symbols}) {
        my $ngene = ($#{$sdat->{hits}} + 1) / 2;
        $counts[$ngene < 10 ? $ngene : 10 ]++;
    }

    my $offSc    = $scH->{uc("Official")};
    my $offAndUn = 0;
    my $multiOff = 0;
    for my $i (0..$#{$ids}) {
        my $dat  = $symbols->{$ids->[$i]};
        my $hits = $dat->{hits};
        ## Tally number of official genes and unofficial genes
        $dat->{o} = 0;
        $dat->{u} = 0;
        for (my $j = 1; $j <= $#{$hits}; $j += 2) {
            my $sc = $hits->[$j];
            $statuses{ $lvls->[ $sc - 1 ] }++;
            if ($sc == $offSc) {
                $dat->{o}++;
            } else {
                $dat->{u}++;
            }
        }
        if (my $oNum = $dat->{o}) {
            $offAndUn++ if ($dat->{u});
            $multiOff++ if ($oNum > 1);
        }
        ## Tally odd characters. Mostly curiosity, but these may
        ## cause issues in some workflows.
        ## Allowed atypical characters:
        ##   '_' eg C4B_2 -> https://www.ncbi.nlm.nih.gov/gene/100293534
        ##   '@' eg HOXA@ -> https://www.ncbi.nlm.nih.gov/gene/3197
        my $sym = $dat->{name};
        $sym =~ s/[a-z0-9\.\-_@]//gi;
        foreach my $char (split('', $sym)) {
            if ($char) {
                $oddChar{$char}{n}++;
                $oddChar{$char}{ex} ||= $dat->{name};
            }
        }
    }

    my $comText = &_factor_map_block( $lvls, $scH, "Nomenclature Status", 
                                      undef, \%statuses );
    $comText .= &_mtx_comment_block("", "'UnofficialPreferred' is simply the first Unofficial symbol listed when no Official symbols are available. It holds no special significance and is provided to allow for one symbol to be recovered consistently for loci that lack an Official symbol.");
    # Summary stats:
    if ($offAndUn || $multiOff) {
        $comText .= &_mtx_comment_block( "Symbols with multiple statuses:","");
        $comText .= sprintf("%% %20s : %6d (Both Official and something else)\n",
           "Official + Not", $offAndUn) if ($offAndUn);
        $comText .= sprintf("%% %20s : %6d ('Official' for 2+ genes! BAD!)\n",
           "Multiple Official", $multiOff) if ($multiOff);
    }

    $comText .= &_mtx_comment_block("", $bar, "",
"CAUTION: Some symbols have case/capitalization subtleties. While there are often historical reasons for these case decisions, attempting to extract information from a gene symbol's case is about as reliable as determining the contents of a wrapped present by listening to the noises it makes when shaken. If you are pivoting data from symbols (rather than accessions) you're already working at a significant disadvantage:",
">PLEASE MATCH SYMBOLS CASE-INSENSITIVELY<",
"The original case is being preserved to aid in 'pretty' display of the symbols.");

    $comText .= "%
% Number of genes refererenced by symbol := Number of symbols
";
    for my $i (1..$#counts) {
       $comText .= sprintf("%%  %s := %d\n", $i == 10 ? ">9" : " $i", $counts[$i] || 0);
    }
    my @ocs = sort { $oddChar{$b}{n} <=> $oddChar{$a}{n} } keys %oddChar;
    if ($#ocs != -1) {
        $comText .= "%
% Potentially troublesome character := Number of symbols (Example)
";
        foreach my $oc (@ocs) {
            my $ocd = $oddChar{$oc};
           $comText .= sprintf("%%  '%s' := %5d (%s)\n", $oc, $ocd->{n}, $ocd->{ex});
        }
    }
    return $comText;
}

sub gene_symbol_meta {
    my ($symbols, $ids, $rc) = @_;
    ## Presumes that &gene_symbol_stats_block( $symbols) has been run
    my @meta = qw(Official NotOfficial);
    my $comText = sprintf("%% %s %s\n", $rc, join($mtxSep, "Name", @meta));
    for my $i (0..$#{$ids}) {
        my $dat  = $symbols->{$ids->[$i]};
        my $hits = $dat->{hits};
        $comText .= sprintf( "%% %d %s\n", $i+1, join($mtxSep, map {
            $dat->{$_} || 0 } qw(name o u)));
    }
    return $comText;
}

sub _generic_meta_block {
    my ($ids, $rc, $meta, $mcols) = @_;
    if (!$meta || !$mcols) {
        $mcols = [];
        $meta  = {};
    }
    my $comText = sprintf( "%% %s %s\n", $rc, join($mtxSep, "Name", @{$mcols}));
    for my $i (0..$#{$ids}) {
        my $id = $ids->[$i];
        my $m = $meta->{$id} || {};

        my @line = ($id, map { defined $_ ? $_ : "" } 
                    map { $m->{$_} } @{$mcols});
        $comText .= sprintf("%% %d %s\n", $i+1, join($mtxSep, @line));
    }
    return $comText;
}

sub _factorize_levels {
    ## Used for GO evidence codes (eclevels) and RefSeq Status (rslevels)
    my $argKey = shift || "";
    my $lvTxt  = $args->{$argKey} || "";
    $lvTxt     = uc($lvTxt) if (shift);
    &death("No levels have been defined for -$argKey") unless ($lvTxt);
    $lvTxt    =~ s/^\s+//; $lvTxt =~ s/\s+$//; # leading/trailing whitespace
    my @lvls  = split(/\s*,\s*/, $lvTxt);
    my %scH  = map { uc($lvls[$_]) => $_ + 1 } (0..$#lvls);
    return (\@lvls, \%scH);
}

## Generates the MTX comment block that encodes factor levels
sub _factor_map_block {
    my ($factorNames, $valueMap, $what, $comments, $counts) = @_;
    my $rv = "%
% $bar
%
%% Values should be treated as factors representing $what:
";

    my $cnt = "%% Count:     ";
    my $top = "%% Index:     ";
    my $bot = "%% LEVELS [,][";
    for my $li (0..$#{$factorNames}) {
        my $l = $factorNames->[$li];
        my $v = $valueMap->{uc($l)};
        my $fmt = '%'.CORE::length($l).'s';
        unless ($li == $#{$factorNames}) {
            $l   .= ',';
            $fmt .= ' ';
        }
        $bot .= $l;
        $top .= sprintf($fmt, $v);
        $cnt .= sprintf($fmt, $counts ? 
                        $counts->{$factorNames->[$li]} || "" : "");
    }
    $rv .= "$cnt\n" if ($counts);
    $rv .= "$top\n";
    $rv .= "$bot]\n";
    map { $rv .= "%% $_\n" } @{$comments} if ($comments);
    return $rv;
}

sub _filter_block {
    my $filt = shift;
    my $rv = &_mtx_comment_block("", $bar, "", "Default Automatic Filters",
 "These parameters define filters that will be automatically applied when the matrix is loaded, unless you set autofilter=FALSE. You can also undo them by calling \$reset() after loading", "");
    
    foreach my $key (sort keys %{$filt}) {
        $rv .= &_default_parameter( $key, $filt->{$key} );
    }
    return $rv;
}

sub _default_parameter {
    my ($key, $val, $com, $sep) = @_;
    return "" unless ($key);
    return "" if (!$val && $val ne "0");
    if (my $r = ref($val)) {
        if ($r eq 'ARRAY') {
            if ($#{$val} == 0) {
                $val = $val->[0];
            } else {
                $sep ||= ',';
                $val = sprintf("[%s][%s]", $sep, join($sep, @{$val}));
            }
        } else {
            die "Unrecognized filter value '$val'";
        }
    }
    my $rv = "";
    $rv .= "% $com\n" if ($com);
    $rv .= sprintf("%%%% DEFAULT %s %s\n", $key, $val);
    return $rv;
}

sub _species_MTX {
    return &_default_parameter("Species", shift);
}

sub _datestamp_for_file {
    my $file = shift;
    if ($file =~ /v(\d{4}-\d{2}-d{2})/) {
        # The file already has a versioned YYYY-MM-DD stamp on it. Use it
        return $1;
    }
    # https://stackoverflow.com/a/1841160
    return POSIX::strftime( "%Y-%m-%d", localtime(  ( stat $file )[9] ) );
}

sub gzfh {
    my ($urlDir, $expectCols) = @_;
    my $src = &fetch_url($urlDir);
    my $fh;
    if ($src =~ /\.gz$/) {
        $fh = IO::Uncompress::Gunzip->new( $src ) ||
            &death("Failed to gunzip file", $src);
    } else {
        open($fh, "<$src") || &death("Failed to read file", $src, $!);
    }
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

## See README.md in the matrixGenerators folder for more information
## on file formats and directory structure

sub file_name {
    my ($type, $mod, $ns1, $ns2, $auth, $vers) = @_;
    my $file = $type.'@';
    $file .= "$mod-" if ($mod);
    $file .= sprintf("%s_to_%s@%s@%s.mtx", $ns1, $ns2, $auth, $vers);
    return $file;
}

sub primary_path {
    return &primary_folder(@_).'/'.&file_name(@_);
}

sub primary_folder {
    my ($auth, $vers);
    if ($#_ == 1) {
        # Just two arguments
        ($auth, $vers) = @_;
    } elsif ($#_ == 5) {
        # Full arguments ($type, $mod, $ns1, $ns2, $auth, $vers)
        ($auth, $vers) = ($_[4], $_[5]);
    }
    my $aDir = sprintf('%s/byAuthority', $outDir );
    my $dir  = sprintf('%s/%s/%s',$aDir, $auth, $vers);
    unless ($tasks{PrimaryFolder}++) {
        ## Make sure folder exists, add README
        &mkpath([$dir]);
        &copy_template_file("byAuthorityReadme.md", "$aDir/README.md");
    }
    return $dir;
}

sub symlinked_paths {
    my ($type, $mod, $ns1, $ns2, $auth, $vers) = @_;
    ## Name of symlink is just the file:
    my $lnk   = &file_name(@_);
    ## Target is up three levels, then into byAuthority:
    my $targ  = sprintf("../../../byAuthority/%s/%s/%s", $auth, $vers, $lnk);
    my $nsDir = sprintf('%s/byNamespace', $outDir );
    unless ($tasks{Symlinks}++) {
        # Set up the namespace folder
         &mkpath([$nsDir]);
         &copy_template_file("byNamespaceReadme.md", "$nsDir/README.md");
    }

    ## Set up namespace hierarchy
    foreach my $np ([$ns1, $ns2], [$ns2, $ns1]) {
        my ($nsa, $nsb) = @{$np};
        my $nd = "$nsDir/$nsa";
        unless ($tasks{"NS-$nsa"}++) {
            &mkpath([$nd]);
            my $rdm = "$nd/README.md";
            unless (-s $rdm) {
                if (open(RDM, ">$rdm")) {
                    print RDM "# Namespace $nsa

This folder holds matrix files that have the `$nsa` namespace assigned
to _either_ their row _or_ column dimensions. The files themselves
will be in one of the subfolders, which describe the other namespace.

";
                    close RDM;
                } else {
                    &err("Failed to create namespace readme", $rdm, $!);
                }
            }
        }

        ## Make the second namespace layer
        $nd = "$nd/$nsb";
        unless ($tasks{"NS-$nsa-$nsb"}++) {
            &mkpath([$nd]);
            my $rdm = "$nd/README.md";
            unless (-s $rdm) {
                if (open(RDM, ">$rdm")) {
                    print RDM "# Namespace Pair $nsa + $nsb

This folder holds matrix files with both the `$nsa` and `$nsb`
namespaces, irrespective of whether one is in rows or columns. The
files are actually symbolic links back to the [byAuthority][BA]
directory.

\[BA\]: ../../../byAuthority

";
                    close RDM;
                } else {
                    &err("Failed to create namespace readme", $rdm, $!);
                }
            }
        }

        ## Make link
        my $nLnk = "$nd/$lnk";
        my $ok = symlink($targ, $nLnk);
        ## We also want to link to the RDS serializations
        symlink("$targ.rds", "$nLnk.rds");
    }

    ### TODO : Find older files, move them to an "archive" folder
}

sub copy_template_file {
    my ($template, $dest, $header) = @_;
    my $tf = "$codeDir/$template";
    
    unless (-s $tf) {
        &err( "Can not copy template '$template' - file not found", $tf);
        return undef;
    }
    if (-s $dest && ((-M $dest) < (-M $tf))) {
        ## Destination file exists and is more recent than template
        return 0;
    }
    if (open(IN, "<$tf")) {
        if (open(OUT, ">$dest")) {
            my $ok = 0;
            while (<IN>) {
                ## Horizontal rule marks start of actual template
                if (/^---/ && !$ok) {
                    $ok = 1;
                    print OUT "$header\n\n" if ($header);
                    next;
                }
                print OUT $_ if ($ok);
            }
            close OUT;
        } else {
            &err("Failed to copy tempate", $dest, $!);
        }
        close IN;
    } else {
        &err("Failed to read template file", $tf, $!);
    }
}

sub parse_filename {
    my $file = shift || "";
    my %rv = (file => $file );
    if ($file =~ /(.+)\.rds$/i) {
        $file = $1;
        $rv{rds} = 1;
    }
    $file =~ s/\.mtx$//i;
    my @bits = split('@', $file);
    if ($#bits == 3) {
        $rv{type} = $bits[0];
        $rv{auth} = $bits[2];
        $rv{vers} = $bits[3];
        my $nstxt = $bits[1];
        if ($nstxt =~ /^(.+)\-(.+)$/) {
            # Modifier / Adjective, eg "Human-Symbol_to_UniProt"
            $rv{mod} = $1;
            $nstxt   = $2;
        }
        my @ns = split('_to_', $nstxt);
        if ($#ns == 1) {
            $rv{rows} = $ns[0];
            $rv{cols} = $ns[1];
        } else {
            $rv{error} = "Namespace '$nstxt' parse failure";
        }
    } else {
        $rv{error} = "Expected 4 fields, found ".($#bits+1);
    }
    return \%rv
}

sub post_process {
    ## Sets symlinks and runs stash, if available
    my $meta = pop @_;
    my $trg  = &primary_path(@_);
    &stash($trg, $meta);
    &symlinked_paths( @_ );
}

sub stash {
    return unless ($stash);
    my $exe;
    if (-s $stash) {
        $exe = $stash;
    } else {
        $exe = `which stash`;
    }
    unless ($exe) {
        warn "
stash is not installed on your system.
    Unless you were really expecting it to be there, this is not an error.

";
        $stash = 0; # Only warn once
        return;
    }
    $exe     =~ s/\s*[\n\r]+$//;
    my ($file, $metaHash) = @_;

    unless ($tasks{StashWarn}++) {
        warn "

Preparing to stash your files...
" ;
        my $sv = `$exe --validate`;
        if ($sv =~ /is valid/) {
            warn "   You are still authenticated\n";
        } else {
            warn "   Please provide your credentials:\n";
            system("$exe --refresh");
        }
        warn "\n";
    }

    unless (-l $file) {
        ## File has not been stashed yet, do so and replace with symlink
        my $cmd = "$exe add --symlink --json \"$file\"";
        my $stashRv = `$cmd`;
    }
    die "Failed to stash file\n  $file\n  " unless (-l $file);
    my $targ = readlink($file);
    my $chksum = basename($targ); $chksum =~ s/\..+$//;

    ## Add/update metadata
    my @meta;
    while (my ($k, $v) = each %{$metaHash || {}}) {
        $k =~ s/[:,]+/_/g;
        my @vals = ($v);
        if (ref($v)) {
            my %u = map { $_ => 1 } @{$v};
            @vals = sort keys %u;
        }
        foreach my $val (@vals) {
            $val =~ s/[:,]+/_/g;
            push @meta, "$k:$val";
        }
    }
    if ($#meta != -1) {
        my $cmd = "$exe metadata --json --metadata '". 
            join(',', @meta)."' $chksum";
        my $stashRv = `$cmd`;
    }
    
    ## Add groups, if requested
    if (my $grp = $args->{stashgroups} || $args->{stashgroup}) {
        my $cmd = "$exe acl --json --groups $grp $chksum";
        my $stashRv = `$cmd`;
    }
    warn "    Stash: $chksum\n";
}


1;

