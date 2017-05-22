#!/usr/bin/perl -w

use strict;
my $scriptDir;
our $defTmp   = "/tmp/pubmedFiles";
my $defFtp   = "ftp.ncbi.nlm.nih.gov";

our $defaultArgs = {
    pubmeddb => "",
    ftp      => $defFtp,
};

BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}
use lib "$scriptDir";
require Utils;
our ($args, $clobber, $ftp, $tmpDir);

use Net::FTP;
use LWP::UserAgent;
use IO::Uncompress::Gunzip;
use XML::Parser::PerlSAX;
use DBI;

# print Dumper($args); die;

my ($sec, $min, $hr, $day, $mon, $year, $wday) = localtime;
my $today = sprintf("%04d-%02d-%02d", $year+1900, $mon+1, $day);

my ($dbh, $xsid, 
    $getXS, $setXS, $doneXS, $clearPM, $delPM, $setPM);
my $defDbName = "simplePubMed.sqlite";
my $maxAbst   = 150;
my $noTitle   = "-No title-";


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

my $fileReq = $args->{file};
my $pmdb    = $args->{pubmeddb};
my $email   = $args->{email};
my $analyze = $args->{analyze};
my $keepxml = $args->{keepxml};

if ($args->{h} || $args->{help} || !($pmdb && $email) ) {
    warn "
Usage:

This program will generate a SQLite database of PubMed IDs, article
titles and publication dates. It was written to support metadata
assignment for PubMed-based gene sets in the SetFisher package. It may
have utility in other applications as well.

The database extracts and stores the following information from all
XML files found in the 'baseline' and 'updatefiles' section of
Medline:

 PubMed ID: The integer PubMed/Medline identifier, eg 17804965

   Version: Almost always 1, but incremented if the same publication
            is present in the system more than once.

      Date: In YYYY-MM-DD format, leaving out day and month if they
            are not identified. Taken from PubDate, or MedlineDate if
            PubDate is not available.

     Title: The article's title, eg 'A Structure for Deoxyribose
            Nucleic Acid'. This will be taken from the XML field
            ArticleTitle, or failing that from VernacularTitle,
            failing that from AbstractText (truncated to $maxAbst
            characters) and finally giving up with '$noTitle'.

XML files will be downloaded from the NCBI FTP site and deleted after
being transformed into TSV files (~10% the size of their XML sources).

Required Arguments:

 -pubmeddb A path to the file that will hold the SQLite
           database. Alternatively, a directory to hold the file, in
           which case the database will be named '$defDbName' at that
           location.

                     >>   ASSURE YOU HAVE SPACE   <<

           In May 2017 the final database was a bit over 5 gigabytes.

    -email NCBI requests that you provide your email address as a
           password when using their FTP servers.

Optional Arguments:

      -ftp Default '$defFtp'. URL to Entrez's FTP site.

   -tmpdir Default '$defTmp'. A location to hold downloaded files

     -file A specific XML file to parse. If not provided, then the
           remote site will be scanned for '*.xml.gz' files, all of
           which will be loaded

   -update If true, indicates that the provided XML file is an update
           (ie, not from the 'baseline' directory). If you forget to
           provide this value, there's a high chance that database
           updates will fail with 'UNIQUE constraint failed' errors.

  -clobber Default 0, which will ignore files already loaded in the
           database. A value of 1 will cause already-loaded XML files
           to be reparsed. A value of 2 will also cause the remote
           file to be fetched again.

  -analyze Default 1, unless -file was specified, in which case
           0. Will ANALYZE the database at the end of the run. Should
           be done periodically to assure statistics are gathered for
           the indices.

  -keepxml Default 0. If true, then XML files will be kept. Otherwise
           they will be deleted after intermediate TSV files have been
           generated from them.

     -help Show this documentation

";
    &err("Please provide an email address for FTP access") unless ($email);
    &err("Please provide the path where you'd like the DB") unless ($pmdb);
    exit;
}

&msg("XML files are: ".($keepxml ?"Kept":"Discarded (set -keepxml to retain)"));

if ($fileReq) {
    &parse_file( $fileReq, $args->{update} );
} else {
    $analyze = 1 unless (defined $analyze);
    &parse_all();
}

if ($analyze) {
    my $t = time;
    &msg("Running ANALYZE on database ... ");
    &get_dbh();
    $dbh->do("ANALYZE");
    warn sprintf("    Done: %.1f min\n", (time - $t) / 60);
}

sub parse_all {
    &msg("Finding XML files at NCBI");
    &get_dbh();
    &_ftp();
    my $tot = 0;
    my @need;
    foreach my $type (qw(baseline updatefiles)) {
        my $sd = "/pubmed/$type";
        $ftp->cwd($sd);
        my @found = map { "$sd/$_" } $ftp->ls('*.xml.gz');
        $tot     += $#found + 1;
        my $upd   = $type eq 'updatefiles' ? 1 : 0;
        
        foreach my $path (@found) {
            my ($xsid, $lastParsed) = &xsid_for_file($path);
            unless ($lastParsed && !$clobber) {
                push @need, [$path, $upd, $xsid, 
                             $lastParsed ? "Reparsing":"New Source"];
            }
        }
    }
    &death("No files found??") unless ($tot);
    my $todo = $#need + 1;
    my $tdf  = $todo / $tot;
    my @bits = (" Total Files: $tot",
                sprintf("Need Parsing: %d = %.1f%%", $todo, 100*$tdf));
                
    if ($tdf > 0) {
        if (my $sz = (-s $pmdb)) {
            push @bits, sprintf("Est. DB Size: %.3f GB", $sz/((1-$tdf)* 10**9));
        }
    }
    &msg(@bits);
    foreach my $fx (@need) {
        &parse_file(@{$fx});
    }
    &msg("All updates finished");
}

sub parse_file {
    my ($file, $isUpd, $xsid, $act) = @_;
    return unless ($file);
    $file =~ s/\.tsv.*$// if ($file =~ /\.tsv/);
    &death("Unexpected file passed for parsing",
           "Files should have format:",
           "medline<stuff>.xml.gz",
           "Request: $file") unless ($file =~ /\bmedline.+\.xml\.gz$/);

    $act ||= "Manual Request";
    unless ($xsid) {
        ## Get the XML Source ID for this file:
        &get_dbh();
        my $lastParsed;
        ($xsid, $lastParsed) = &xsid_for_file($file);
        if ($lastParsed) {
            ## We have already parsed
            if ($clobber) {
                $act = "Reparsing";
            } else {
                &msg("Already parsed: $file");
                return;
            }
        } else {
            $act = "New Source";
        }
    }
    &msg("$act: $file");
    my $tsv = &make_tsv($file); 
    return unless ($tsv);
    my $t = time;
    ## Clear any old records for this XML file:
    $clearPM->execute( $xsid );
    open(TSV, "<$tsv") || &death("Failed to read TSV file", $tsv, $!);
    my $head = <TSV>;
    my $num  = 0;
    if (!defined $isUpd && $file =~ /update/) {
        &msg("Presuming that requested file is an update");
        $isUpd = 1;
    }
    while (<TSV>) {
        s/[\n\r]+$//;
        my ($pmid, $v, $date, $title) = split(/\t/);
        ## Remove any old entries if this is an update file:
        $delPM->execute($pmid, $v) if ($isUpd);
        $setPM->execute($pmid, $v, $xsid, $date, $title);
        $num++;
    }
    close TSV;
    $doneXS->execute($today, $xsid);
    $dbh->commit;
    my @bits = ("Entries: $num (xsid = $xsid)",
                sprintf("DB Time: %d sec", time - $t));
    &msg(@bits);
}

sub make_tsv {
    my ($fReq) = @_;
    my $tsv = join('/', $tmpDir, basename($fReq). ".tsv");
    return $tsv unless (&output_needs_creation($tsv));
    my $file = &fetch_url($fReq);
    unless (-s $file) {
        &err("Failed to recover file, skipping:",$_[0]);
        return "";
    }

    my $t = time;

    my $handler = PubMedHandler->new( $tsv );
    my $parser  = XML::Parser::PerlSAX->new( Handler => $handler );
    my $fh;
    if ($file =~ /\.gz$/) {
        $fh = IO::Uncompress::Gunzip->new( $file ) ||
            &death("Failed to gunzip file", $file);
    } else {
        open($fh, "<$file") || &death("Failed to read file", $file, $!);
    }
    
    # $dbh->begin_work;
    eval {
        $parser->parse( Source => { ByteStream => $fh } );
    };
    if ($handler->{doneXML}) {
        my @bits = (" TSV: $tsv",
                    sprintf("Time: %.1f min", (time - $t) / 60));
        unless ($keepxml) {
            unlink($file);
        }
        &msg(@bits);
    } else {
        unlink($tsv);
        &death("Failed to complete parse of XML", $file);
    }
    return $tsv;
}


sub get_dbh {
    return $dbh if ($dbh);
    if (-d $pmdb) {
        $pmdb =~ s/\/+$//;
        $pmdb = "$pmdb/$defDbName";
    }
    return &_create_db( $pmdb ) unless ( -s $pmdb );
    return &_sths( $pmdb );
 }

sub _con_db {
    my $file = shift;
    &msg("Connecting to DB: $file");
    return $dbh = 
        DBI->connect("dbi:SQLite:dbname=$file",'','', {
            AutoCommit => 0,
            RaiseError => 1,
            PrintError => 0 });
}

sub _sths {
    my $file = shift;
    $dbh ||= &_con_db($file);

    $getXS   = $dbh->prepare("SELECT xsid, parsed FROM xml_source".
                             " WHERE file = ?");
    $setXS   = $dbh->prepare("INSERT INTO xml_source (file) VALUES (?)");
    $doneXS  = $dbh->prepare("UPDATE xml_source SET parsed = ? WHERE xsid = ?");

    $clearPM = $dbh->prepare("DELETE FROM pmid WHERE xsid = ?");
    $delPM   = $dbh->prepare("DELETE FROM pmid WHERE pmid = ? AND vers = ?");
    $setPM   = $dbh->prepare
        ("INSERT INTO pmid (pmid, vers, xsid, pubdate, title)".
         " VALUES (?, ?, ?, ?, ?)");
    return $dbh;
}

sub xsid_for_file {
    my ($file) = @_;
    return wantarray ? (0,"") : 0 unless ($file);
    my $short = basename($file);
    $getXS->execute( $short );
    my $rv = $getXS->fetchall_arrayref();
    if ($#{$rv} == -1) {
        ## Need to note this file in the DB
        eval {
            $setXS->execute( $short);
            $dbh->commit;
        };
        $getXS->execute( $short );
        $rv = $getXS->fetchall_arrayref();
        &death("Failed to insert new XML source", $file) if ($#{$rv} != 0);
    } elsif ($#{$rv} != 0) {
        &death("XML source is not uniquely represented in DB", $file);
    }
    return wantarray ? @{$rv->[0]} : $rv->[0][0];
}


sub _create_db {
    my $file = shift;
    &msg("Creating SQLite database", $file);
    &_con_db($file);
    # print Dumper( $dbh->sqlite_db_status() );
    ## Table to track files being parsed

    ## Primary integer keys in SQLite act as sorta-kinda
    ## auto-incrememnt fields, with the danger of value recycling
    ## after delete:
    ## https://stackoverflow.com/a/7906029

    $dbh->do( "CREATE TABLE xml_source (
  xsid   INTEGER PRIMARY KEY,
  file   TEXT,
  parsed TEXT
)");
    $dbh->do("CREATE UNIQUE INDEX xf_idx ON xml_source (file)");
    
    ## xsid     INTEGER,
    ## FOREIGN KEY(xsid) REFERENCES xml_source(xsid),
    $dbh->do( "CREATE TABLE pmid (
  pmid     INTEGER,
  vers     INTEGER,
  xsid     INTEGER,
  pubdate  TEXT,
  title    TEXT,
  FOREIGN KEY(xsid) REFERENCES xml_source(xsid)
)");
    ## A PMID can exist in multiple versions, and those versions can
    ## be in multiple files. This makes uniqueness management a bit more
    ## awkward; Can't handle it simply file-by-file. Example:
    ## PMID:27928497 in 17n0891 + 17n0892
    $dbh->do( "CREATE UNIQUE INDEX pmid_unique_idx ON pmid (pmid, vers)");
    $dbh->commit;
    
    return &_sths( $file );
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###   Custom SAX handler module
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

package PubMedHandler;

use Data::Dumper;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        stack  => [],
        record => {},
        done   => 0,
    };
    bless ($self, $class);
    my $tsv = $self->{tsv} = shift;
    my $tmp = $self->{tmp} = "$tsv.tmp";
    my $fh;
    open($fh, ">$tmp") || &death("Failed to write TSV", $tmp, $!);
    print $fh join("\t", qw(pmid vers pubdate title))."\n";
    $self->{fh} = $fh;
    return $self;
}

sub DESTROY {
    my $self = shift;
    unless ($self->{doneXML}) {
        &err("Object destruction before XML parse completion!");
        if (my $fh = $self->{fh}) {
            close $fh;
        }
        # Do not allow partial TSV file:
        if (my $tsv = $self->{tsv}) {
            unlink($tsv) if (-s $tsv);
        }
    }
}

sub msg {
    return &main::msg(@_);
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
    if (my $meth = $self->can("_end_$name")) {
        &{$meth}($self, $name, $node);
    }
}

sub _end_PubmedArticle {
    my ($self, $name, $node) = @_;
    # Got to the end of a record, write it
    $self->_write_record();
}

sub _end_PubmedArticleSet {
    my ($self, $name, $node) = @_;
    # Full document was parsed
    $self->{doneXML} = 1;
    if (my $fh = $self->{fh}) {
        close $fh;
        my $tsv = $self->{tsv};
        my $tmp = $self->{tmp};
        rename($tmp, $tsv);
    }
}

sub _end_PMID {
    my ($self, $name, $node) = @_;
    # die Dumper($node);
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
            $self->{record}{idv}  = -2;
            if (my $att = $node->{Attributes}) {
                $self->{record}{idv} = $att->{Version} || -1;
            }
        }
    }
}

# Globs do not seem to work with ->can("globname")
sub _end_ArticleTitle    { shift->_end_Title_set( @_ ) }
sub _end_VernacularTitle { shift->_end_Title_set( @_ ) }
sub _end_Title_set {
    my ($self, $name, $node) = @_;
    # Hopefully the title, which is 99% of why I am interested in
    # parsing these files. Some entries have empty ArticleTitle, but
    # filled VernacularTitle, eg PMID:27868948
    if (my $txt = $node->{text}) {
        ## These tags can be empty (no text) eg PMID:27868948
        push @{$self->{record}{article}[-1]{$name}}, join('', @{$txt});
    }    
}

sub _end_AbstractText {
    my ($self, $name, $node) = @_;
    ## Using this as a fallback title text
    if (my $txt = $node->{text}) {
        $txt = join('', @{$txt});
        ## Truncate to $maxAbst characters
        $txt = substr($txt, 0, $maxAbst).'...' 
            if (CORE::length($txt) > $maxAbst );
        ## Abstract sections can have multiple entries. Just take the
        ## first encounted one, generally will be introduction or
        ## background. eg PMID:27735882
        $self->{record}{article}[-1]{Abstract} ||= [ $txt ];
    }
}
 
sub _end_Year        { shift->_end_Date_set( @_ ) }
sub _end_Month       { shift->_end_Date_set( @_ ) }
sub _end_Day         { shift->_end_Date_set( @_ ) }
sub _end_MedlineDate { shift->_end_Date_set( @_ ) }
sub _end_Date_set {
    my ($self, $name, $node) = @_;
    # Date information (Year, Month, Day and grab-bag "MedlineDate")
    my $par = $self->{stack}[-1]{Name};
    if ($par eq 'PubDate') {
        # Date information shows up many places, make sure it is
        # attached to the article
        push @{$self->{record}{article}[-1]{$name}}, join('', @{$node->{text}});
        #if ($name eq 'MedlineDate') { warn Dumper($self->{record}); die; }
    }
}

sub _end_Title {
    my ($self, $name, $node) = @_;
    my $par = $self->{stack}[-1]{Name};
    if ($par eq 'Journal') {
        # This appears to be the title of the journal. I am not doing
        # anything with these at the moment. If someone were
        # interested in getting this information, be sure to check if
        # other parent tags (eg maybe "Conference", "TextBook",
        # "SketchyBlog", etc) might be present
        push @{$self->{record}{article}[-1]{$par}}, join('', @{$node->{text}});
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
    return unless ($record);
    my $pmid = $record->{pmid};
    return unless ($pmid);
    open(FOO, ">/tmp/foo.txt"); print FOO "Starting PMID $pmid\n"; close FOO;
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
    my $title = $art->{ArticleTitle} || $art->{VernacularTitle} || 
        $art->{Abstract};
    unless ($title) {
        # No title found
        if (my $jour = $art->{Journal}) {
            # At least note the journal
            $title = [ map { "[Untitled] $_" } @{$jour} ] if ($jour->[0]);
            ## ... but these seem to get updated! eg PMID:27869070,
            ## which is untitled in 2016 XML but has title ("Une cause
            ## inhabituelle d'hyperCKemie.") at NLM as of May 2017.
        } else {
            &msg("No title found for PMID:$pmid");
        }
        $title ||= [$noTitle];
    }
    &msg("Multiple titles for PMID:$pmid") if ($#{$title} > 0);
    $title = $title->[0] || "";
    ## Life is much easier if we stick to ASCII. This will remove some
    ## special symbols
    $title =~ s/\P{IsASCII}/?/g;
    $title =~ s/\s*\t\s*/ /g; # Tabs would be A Bad Thing
    my $fh = $self->{fh};
    print $fh join
        ("\t", $pmid, $record->{idv} || 0, $date, $title)."\n";
    unless (++$self->{done} % 5000) {
        &msg(sprintf("  #%5d %15s [%-10s] %s", $self->{done}, "PMID:$pmid",
                     $date, substr($title, 0, 80) ) );
    }
    # printf(" %8s [%s] %s\n", "PMID:$pmid", $date, $title->[0]);
}
