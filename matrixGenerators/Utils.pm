use strict;
use Data::Dumper;
use File::Path 'mkpath';
use Net::FTP;
use LWP::UserAgent;

our ($defaultArgs, $defTmp, $ftp, $ua);

our $args = &parseargs({
    ## Tool-specific defaults:
    %{$defaultArgs},
    ## Defaults common to tools:
    tmpdir   => $defTmp,
    clobber  => 0,
    verify_hostname => 0,  });

our $tmpDir   = $args->{tmpdir}; $tmpDir =~ s/\/+$//;
our $clobber  = $args->{clobber} || 0;
&mkpath([$tmpDir]) if ($tmpDir);

# die Dumper($args);

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

sub _ftp {
    ## http://perlmeme.org/faqs/www/ftp_file_list.html
    return $ftp if ($ftp);
    my $site = $args->{ftp};
    $ftp = Net::FTP->new($site) ||
        &death("Failed to connect to FTP", $site, $!);
    my $pass = $site =~ /ncbi/ ? $args->{email} : "";
    $ftp->login("anonymous", $pass) ||
        &death("Failed to login to FTP", $site, $!);
    $ftp->binary();
    return $ftp;
}

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

sub fetch_url {
    my ($uReq, $dReq) = @_;
    my $dest = &local_file( @_ );
    if (&source_needs_recovery($dest)) {
        ## File not yet downloaded, or request to re-download
        my $pending = 1;
        ## For NCBI it looks like the timeout is 60 seconds, and
        ## parsing is ~90 seconds per file on my system. Close and
        ## reinitialize to be assured of a 'live' connection:
        undef $ftp; &_ftp();
        while ($pending) {
            $ftp->get($uReq, $dest);
            if (-s $dest) {
                &msg("Downloaded $uReq:", $dest);
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
    my ($url, $dReq) = @_;
    my $dest = &local_file( @_ );
    if (&source_needs_recovery($dest)) {
        &_ua();
        my $resp = $ua->get($url, ':content_file' => $dest );
        if ($resp->is_success) {
            if (-s $dest) {
                &msg("Downloaded $url:", $dest);
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


1;
