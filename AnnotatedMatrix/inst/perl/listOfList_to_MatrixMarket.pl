#! /usr/bin/perl -w

use strict;

=pod

Converts an internal "list of list" formats into a MatrixMarket file
for ease of reading in R.

=cut
    
my $file = $ARGV[0];

die "Please provide the file to convert as the first argument\n" unless ($file);

my $out = $ARGV[1] || "$file.mtx";

warn "Input file: $file\n";

open(IN, "<$file") || die "Failed to read input\n  $file\n  $!";
open(OUT,">$out")  || die "Failed to write output\n  $out\n  $!";
my (%listHash, %listCount, %rowOrder, %colOrder);
my ($colCount, $rowCount) = (0,0);

my $lname = "";
my $tot   = 0;

my $listsAsCols = 1;

my $setFunc = $listsAsCols ? sub {
    my ($lname, $id) = @_;
    my $rid = $rowOrder{$id}    ||= ++$rowCount;
    my $cid = $colOrder{$lname} ||= ++$colCount;
    unless ($listHash{$rid}{$cid}) {
        $listHash{$rid}{$cid} = ++$listCount{$lname};
        $tot++;
    }
} : sub {
    my ($lname, $id) = @_;
    my $rid = $rowOrder{$lname} ||= ++$rowCount;
    my $cid = $colOrder{$id}    ||= ++$colCount;
    unless ($listHash{$rid}{$cid}) {
        $listHash{$rid}{$cid} = ++$listCount{$lname};
        $tot++;
    }
};

my (%kv, %metadata);

while (<IN>) {
    s/[\n\r]+$//;
    if (/^#/) {
        # Comment row
        if (/^#+\s+LIST\s+-\s+(.+?)\s*$/) {
            $lname = $1;
        } elsif (/^\#+\s*(\S+?)=\s*(.+?)\s*$/) {
                 # keyval pair
            if ($lname) {
                $metadata{$1}{$lname} = $2;
            } else {
                $kv{$1} = $2;
            }
        }
        next;
    }
    next if (/^\s*$/); # Blank line
    die "Identifier ($_) encountered before list name was set\n  "
        unless $lname;
    &{$setFunc}( $lname, $_ )
}
close IN;

my @rows = sort { $rowOrder{$a} <=> $rowOrder{$b} } keys %rowOrder;
my @cols = sort { $colOrder{$a} <=> $colOrder{$b} } keys %colOrder;

my $bar = "%-------------------------------------------\n";
my $sep = " :: ";
print OUT "%%MatrixMarket matrix coordinate integer general\n";
print OUT "% Separator '$sep'\n";

foreach my $k (sort keys %kv) {
    printf(OUT "%% DEFAULT %s %s\n", $k, $kv{$k});
}

print OUT $bar;
my @mc = sort keys %metadata;
my @rhead = ("Row Name");
my @chead = ("Col Name");
if ($listsAsCols) {
    push @chead, @mc;
} else {
    push @rhead, @mc;
}
printf(OUT "%% %s\n", join($sep, @rhead));
for my $i (0..$#rows) {
    my $id = $rows[$i];
    my @r = ($id);
    if (!$listsAsCols) {
        push @r, map { $metadata{$_}{$id} || "" } @mc;
    }
    printf(OUT "%% %d %s\n", $i + 1, join($sep, @r));
}

print OUT $bar;
printf(OUT "%% %s\n", join($sep, @chead));
for my $i (0..$#cols) {
    my $id = $cols[$i];
    my @r = ($id);
    if ($listsAsCols) {
        push @r, map { $metadata{$_}{$id} || "" } @mc;
    }
    printf(OUT "%% %d %s\n", $i + 1, join($sep, @r));
}

print OUT $bar;
print OUT "% Matrix triples : List Identifier Order\n";
print OUT $bar;
printf(OUT "  %d %d %d\n", $#rows + 1, $#cols + 1, $tot);
foreach my $lid (sort { $a <=> $b } keys %listHash) {
    my $iH = $listHash{$lid};
    foreach my $iid (sort { $a <=> $b } keys %{$iH}) {
        printf(OUT "%d %d %d\n", $lid, $iid, $iH->{$iid});
    }
}

close OUT;
warn "Output file: $out\n";
