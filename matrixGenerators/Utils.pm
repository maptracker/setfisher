use strict;
use Data::Dumper;
use File::Path 'mkpath';

our ($defaultArgs, $defTmp, $ftp);

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



sub fetch_url {
    my ($uReq, $dReq) = @_;
    unless ($dReq) {
        $dReq = $uReq; $dReq =~ s/.+\///;
    }
    my $dest = "$tmpDir/$dReq"; # Local file path
    if (&source_needs_recovery($dest)) {
        ## File not yet downloaded, or request to re-download
        my $pending = 1;
        while ($pending) {
            $ftp->get($uReq, $dest);
            if (-s $dest) {
                &msg("Downloaded $dReq", $dest);
                $pending = 0;
            } else {
                if ($pending >= 3) {
                    &err("Failed to recover remote file after $pending tries",
                         "Source: $uReq", "Destination: $dest", $@);
                    return "";
                }
                $pending++;
            }
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
