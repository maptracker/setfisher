#!/usr/bin/perl -w

use strict;
use LWP::UserAgent;
use IO::Uncompress::Gunzip;
use XML::Parser::PerlSAX;

use Data::Dumper;

my @mnThree    = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @mnFull     = qw(January February March April May June 
                    July August September October November December);
# These are not quite right, but close enough:
my $mnHash     = {
    winter    => 1,
    spring    => 4,
    spr       => 4,
    summer    => 7,
    sum       => 7,
    autumn    => 10,
    fall      => 10,
    easter    => 4,  # http://www.ncbi.nlm.nih.gov/pubmed/2019856
    christmas => 12, # http://www.ncbi.nlm.nih.gov/pubmed/13954194
    '1st Quart' => 1, # PMID:10237212
    '2d Quart'  => 4, # PMID:10248150
    '3d Quart'  => 7, # PMID:10236507
    '4th Quart' => 10, # PMID:10249456
    '-00 Winter' => 1, # PMID:10711319
    '-94 Winter' => 1, # PMID:11362190
    'N0v' => 11, # PMID:5275180
    '' => 1, #
    '' => 1, #
    '' => 1, #
    '' => 1, #
    '' => 1, #
    '' => 1, #
};
foreach my $arr (\@mnThree, \@mnFull) {
    map { $mnHash->{lc($arr->[$_])} = $_ + 1 } (0..$#{$arr});
}
for my $i (1..12) {
    $mnHash->{$i} = $i;
    $mnHash->{sprintf("%02d", $i)} = $i;
}
my $grabDate = { map { $_ => 1 } qw(Year Month Day MedlineDate) };



&parse_file($ARGV[0]);

sub parse_file {
    my $file = shift;
    return unless ($file);

    my $handler = PubMedHandler->new(  );
    my $parser  = XML::Parser::PerlSAX->new( Handler => $handler );
    my $fh;
    if ($file =~ /\.gz$/) {
        $fh = IO::Uncompress::Gunzip->new( $file ) ||
            &death("Failed to gunzip file", $file);
    } else {
        open($fh, "<$file") || &death("Failed to read file", $file, $!);
    }
    
    $parser->parse( Source => { ByteStream => $fh } );

}

sub msg {
    warn "[*] ".join("\n    ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}

sub err {
    warn "[!!] ERROR: ".join
        ("\n     ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}

sub death { &err(@_); die " -- "; }

### Custom SAX handler module:

package PubMedHandler;

use Data::Dumper;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        stack => [],
        record => {},
    };
    bless ($self, $class);
    return $self;
}

sub start_element {
    my ($self, $element) = @_;
    my $name = $element->{Name};
    push @{$self->{stack}}, $element;
    if ($name eq 'PubmedArticle') {
        $self->{record} = {};
    } elsif ($name eq 'Article') {
        push @{$self->{record}{article}}, {};
    }

}

sub end_element {
    my ($self, $element) = @_;
    my $node = pop @{$self->{stack}};
    my $name = $node->{Name};
    if ($name eq 'PubmedArticle') {
        # Got to the end of a record, write it
        $self->_write_record();
    } elsif ($name eq 'PMID') {
        # Primary ID
        my $par = $self->{stack}[-1]{Name};
        if ($par eq 'MedlineCitation') {
            # PMID entries can show up elsewhere, so check parent
            my $pmid = join('', @{$node->{text}});
            ## die "DEBUGGING" if ($pmid > 100);
            if ($self->{record}{pmid}) {
                &err("Multiple PMID entries", "$pmid vs $self->{record}{pmid}");
            } else {
                $self->{record}{pmid} = $pmid;
            }
        }
    } elsif ($name eq 'ArticleTitle') {
        # Hopefully the title, which is 99% of why I am interested in
        # parsing these files.
        push @{$self->{record}{article}[-1]{Title}}, 
        join('', @{$node->{text}});
    } elsif ($grabDate->{$name}) {
        # Date information (Year, Month, Day and grab-bag "MedlineDate")
        my $par = $self->{stack}[-1]{Name};
        if ($par eq 'PubDate') {
            # Date information shows up many places, make sure it is
            # attached to the article
            push @{$self->{record}{article}[-1]{$name}}, 
            join('', @{$node->{text}});
            #if ($name eq 'MedlineDate') { warn Dumper($self->{record}); die; }
        }
    } elsif ($name eq 'Title') {
        my $par = $self->{stack}[-1]{Name};
        if ($par eq 'Journal') {
            # This appears to be the title of the journal. I am not
            # doing anything with these at the moment. If someone were
            # interested in getting this information, be sure to check
            # if other parent tags (eg maybe "Conference", "TextBook",
            # "SketchyBlog", etc) might be present
            push @{$self->{record}{article}[-1]{$par}}, 
            join('', @{$node->{text}});
        }
    }
}

sub characters {
    my ($self, $chars) = @_;
    push @{$self->{stack}[-1]{text}}, $chars->{Data};
    # Do nothing if we are not actively parsing a target
    #return unless ($self->{ACTIVE});
    #my $parent = $self->{STACK}[-1];
    #my $txt    = $chars->{Data};
    #push @{$parent->{TEXT}}, $txt;
}

sub _write_record {
    my $self = shift;
    my $record = $self->{record};
    #print Dumper($record); die;
    my $pmid = $record->{pmid};
    return unless ($record);
    my $date = "";
    my @arts = @{$record->{article} || []};
    &msg("Multiple articles for PMID:$pmid") if ($#arts > 0);
    my $art = $arts[0] || {};
    if (my $y = $art->{Year}) {
        &msg("Multiple years for PMID:$pmid") if ($#{$y} > 0);
        $y = $y->[0];
        if ($y =~ /^\d{4}$/) {
            $date = $y;
            if (my $m = $art->{Month}) {
                &msg("Multiple months for PMID:$pmid") if ($#{$m} > 0);
                $m = $m->[0];
                if (my $std = $mnHash->{lc($m)}) {
                    $date = sprintf("%s-%02d", $date, $std);
                    if (my $d = $art->{Day}) {
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
    } elsif (my $mld = $art->{MedlineDate}) {
        ## eg: <MedlineDate>1978 Sep-Oct</MedlineDate> This appears to
        ## be a fallback value, AFAICT used when a date spans a range
        ## rather than having a discrete value.
        &msg("Multiple Medline dates for PMID:$pmid") if ($#{$mld} > 0);
        if ($mld->[0] =~ /^(\d{4})\b/) {
            ## Do not trust getting more than a year here
            $date = $1;
        }
    }
    my $title = $art->{Title} ||= ["-No title-"];
    &msg("Multiple titles for PMID:$pmid") if ($#{$title} > 0);
    ## Life is much easier if we stick to ASCII
    $title->[0] =~ s/\P{IsASCII}/?/g;
    printf(" %8s [%s] %s\n", "PMID:$pmid", $date, $title->[0]);
}

sub msg {
    warn "[*] ".join("\n    ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}

sub err {
    warn "[!!] ERROR: ".join
        ("\n     ", map { defined $_ ? $_ : '-UNDEF-' } @_). "\n";
}
