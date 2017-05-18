use strict;
use LWP::UserAgent;
use Data::Dumper;

our ($defaultArgs, $defTmp, $ftp);

our $args = &parseargs({
    ## Tool-specific defaults:
    %{$defaultArgs},
    ## Defaults common to tools:
    tmpdir   => $defTmp,
    clobber  => 0,
    verify_hostname => 0,
                       });

my $tmpDir   = $args->{tmpdir}; $tmpDir =~ s/\/+$//;
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



sub fetch_url {
    my ($uReq, $dReq) = @_;
    unless ($dReq) {
        $dReq = $uReq; $dReq =~ s/.+\///;
    }
    my $dest = "$tmpDir/$dReq"; # Local file path
    if (&source_needs_recovery($dest)) {
        ## File not yet downloaded, or request to re-download
        my $url = "$ftp/$uReq"; # Remote URL
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


1;
