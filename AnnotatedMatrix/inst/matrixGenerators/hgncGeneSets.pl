#!/usr/bin/perl -w


use strict;
my $scriptDir;
our $defTmp = "/tmp/hgncGeneSets";
my $defFtp = "ftp.ebi.ac.uk";
my $defSym = "Unknown,Unofficial,UnofficialPreferred,Interim,Official";
my $gnDom  = "http://www.genenames.org";

my $symUrl    = $gnDom.'/cgi-bin/gene_search?search=%s';
my $hgncUrl   = $gnDom.'/cgi-bin/gene_symbol_report?hgnc_id=%s';
my $entrezUrl = 'https://www.ncbi.nlm.nih.gov/gene/%s'; # For integer IDs
my $omimUrl   = 'https://omim.org/entry/%s';
my $ensgUrl   = 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s';

my %nsUrl = ( HGNC    => $hgncUrl,
              Entrez  => $entrezUrl,
              Ensembl => $ensgUrl,
              OMIM    => $omimUrl);

=head Random Things

=head2 Column: status

Values Approved or Entry Withdrawn

=head2 Column: gene_family

~2000 gene categories. Most (60%)have just one gene. Maybe useful as ontology?

=head2 Column: locus_group

  19098 protein-coding gene
  13077 pseudogene
   7145 non-coding RNA
   1205 other
   1161 withdrawn
    578 phenotype

Column locus_type is more-or-less the same, but a bit more detailed
(29 categories)

=cut

our $defaultArgs = {
    ftp       => $defFtp,
    symlevels => $defSym,
    dir       => ".",
};

BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}
use lib "$scriptDir";
require Utils;
our ($args, $clobber, $ftp, $tmpDir, $maxAbst, $bar);

my $outDir   = $args->{dir};    $outDir =~ s/\/+$//;

&mkpath([$outDir]);


if ($args->{h} || $args->{help}) {
    warn "
Usage:

This program will generate gene sets for the SetFisher enrichment
package. It will recover information from the HUGO Gene Nomenclature
Committee (HGNC) FTP site and build annotated MatrixMarket files from
the data recovered there.

Optional Arguments:

      -dir Default '.'.  Output directory that will contain generated
           matrix and metadata files

  -clobber Default 0, which will preserve any already-generated
           files. Set to 1 to regenerate matrix files, and any value
           greater than 1 to also re-download base files from Entrez.

";
    exit;
}

my $src  = "pub/databases/genenames/new/tsv/hgnc_complete_set.txt";
my $loc  = &fetch_url($src);
my $vers = &_datestamp_for_file($loc);

warn "

Building pivot matrices from HGNC:

   FTP Site: $defFtp
     Source: $src
 Local Copy: $loc
  Datestamp: $vers

";

&symbol_maps();

sub symbol_maps {
    # Grab the following namespaces:
    my ($rows, $xcol) =
        &_get_rows([ HGNC    => 'hgnc_id', # Easier to grab this value twice
                     Entrez  => 'entrez_id',
                     Ensembl => 'ensembl_gene_id' ]);
    
    ## Factorize the symbols
    my ($lvls, $scH) = &_factorize_levels('symlevels');
    my $offL = $scH->{ uc("Official") };
    my $upL  = $scH->{ uc("UnofficialPreferred") };
    my $uL   = $scH->{ uc("Unofficial") };
    my $oldL = $scH->{ uc("Unknown") };
    foreach my $row (@{$rows}) {
        my %syms;
        if (my $off = $row->[0][1]) {$syms{$off} ||= $offL;} # Official
        if (my $pre = $row->[0][4]) {$syms{$pre} ||= $upL;}  # UnPreferred
        map { $syms{$_} ||= $uL   } @{$row->[0][5]}; # All other Unofficial
        map { $syms{$_} ||= $oldL } @{$row->[0][6]}; # "Previous" -> Unknown
        $row->[2] = \%syms;
    }

    ## Now build the lookups for the 'main' namespaces
    my @gMetaCols = qw(Symbol Type Description);
    for my $xi (0..$#{$xcol}) {
        my $ns  = $xcol->[$xi];
        my $trg = sprintf("%s/Map-HGNC_Symbol-to-%s-v%s.mtx", 
                      $outDir, $ns, $vers);
        unless (&output_needs_creation($trg)) {
            ## File is already made, and clobber is not set
            &msg("Keeping existing Symbol map:", $trg);
            next;
        }

        my (%symbols, %genes, %geneMeta);
        my ($gnum, $nznum) = (0, 0);
        foreach my $row (@{$rows}) {
            my $xid = $row->[1][$xi];
            next unless ($xid);
            $geneMeta{$xid} = {
                Symbol => $row->[0][1] || $row->[0][4] || "",
                Type   => $row->[0][2],
                Description => $row->[0][3],
            };
            my $gn ||= $genes{$xid} ||= ++$gnum;
            while (my ($sym, $sc) = each %{$row->[2]}) {
                ## We are going to preserve case!! ABC != Abc
                my $targ = $symbols{$sym} ||= {
                    name => $sym,
                    hits => [],
                };
                push @{$targ->{hits}}, ($gn, $sc);
                $nznum++;
            }
        }

        my @rowIds  = sort keys %symbols;
        my @gids    = sort { $genes{$a} <=> $genes{$b} } keys %genes;
        my ($rnum,$cnum) = ($#rowIds + 1, $#gids + 1);


        my $tmp = "$trg.tmp";
        open(MTX, ">$tmp") || &death("Failed to write Symbol map", $tmp, $!);
        print MTX &_initial_mtx_block
            ( "Mapping", $rnum, $cnum, $nznum, "HGNC Symbol-to-$ns Map",
              "Scores are factor levels representing the nomenclature status of the assignment.",
              "Gene Symbol", $symUrl,
              "$ns ID", $nsUrl{$ns});
        print MTX &gene_symbol_stats_block( \%symbols, \@rowIds, $lvls, $scH );
        print MTX &_rowcol_meta_comment_block();

        ## Note row metadata
        print MTX &gene_symbol_meta( \%symbols, \@rowIds, 'Row');

        ## Column metadata
        print MTX "% $bar\n";
        print MTX &_generic_meta_block(\@gids, 'Col', \%geneMeta, \@gMetaCols);

        print MTX &_triple_header_block( $rnum, $cnum, $nznum );
        for my $i (0..$#rowIds) {
            my $hits = $symbols{$rowIds[$i]}{hits};
            for (my $j = 0; $j < $#{$hits}; $j += 2) {
                printf(MTX "%d %d %0.1f\n", $i+1, $hits->[$j], $hits->[$j+1]);
            }
        }
        close MTX;

        &msg("Generated Symbol to $ns mapping", $trg);

    }

}

sub _get_rows {
    # Default columns we'll get for all maps:
    my $cols = { 
        hgnc  => 'hgnc_id',
        sym   => 'symbol',
        alias => 'alias_symbol',
        old   => 'prev_symbol',
        desc  => 'name',
        type  => 'locus_group',
    };
    # Additional columns of interest:
    my $xtra = shift;
    my @extra;
    for (my $i = 0; $i < $#{$xtra}; $i += 2) {
        my ($tag, $col) = ($xtra->[$i], $xtra->[$i+1]);
        $cols->{$tag} = $col;
        push @extra, $tag;
    }

    my ($fh) = &gzfh($src, $cols);
    my $hInd = $cols->{hgnc};  # HGNC ID, eg HGNC:15550
    my $sInd = $cols->{sym};   # "Main" symbol, eg ARID4B
    my $aInd = $cols->{alias}; # Alias symbol, eg BCAA|BRCAA1|SAP180
    my $oInd = $cols->{old};   # "Previous" symbol, eg RBP1L1
    my $tInd = $cols->{type};  # Gene type, eg protein-coding gene
    my $dInd = $cols->{desc};  # Description, eg AT-rich interaction domain 4B

    my @xinds = map { $cols->{$_} } @extra;
    
    my $rv = [];
    while (<$fh>) {
        s/[\n\r]+$//;
        my @row = split("\t");
        @row = map { /\"(.+)\"/ || /\'(.+)\'/ ? $1 : $_ } @row;
        ## HGNC, Official symbol, Type, Description
        my @basic = map { $row[ $_ ] } ($hInd, $sInd, $tInd, $dInd);
        ## Aliases, and picking firsts as "UnofficialPreferred" if no Official:
        my $alias = [$row[$aInd] ? split(/\|/, $row[$aInd]) : () ];
        push @basic, $basic[1] ? "" : shift @{$alias} || "";
        push @basic, $alias;
        ## Previous symbols
        push @basic, [$row[$oInd] ? split(/\|/, $row[$oInd]) : () ];
        my @x = map { $row[$_] } @xinds;
        push @{$rv}, [\@basic, \@x];
    }
    return ($rv, \@extra);
}

