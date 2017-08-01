#!/usr/bin/perl -w

use strict;
my $scriptDir;
our $defTmp = "/tmp/hgncGeneSets";
my $defFtp = "ftp.ebi.ac.uk";
my $defSym = "Unknown,Unofficial,UnofficialPreferred,Interim,Official";
my $defTyp = "Deprecated,Unknown,Phenotype,Pseudogene,Other,ncRNA,ProteinCoding";
my $gnDom  = "http://www.genenames.org";

my $entrezUrl = 'https://www.ncbi.nlm.nih.gov/gene/%s'; # For integer IDs
my $omimUrl   = 'https://omim.org/entry/%s';
my $ensgUrl   = 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s';

my $nsUrl = {
    Symbol     => $gnDom.'/cgi-bin/gene_search?search=%s',
    HGNC       => $gnDom.'/cgi-bin/gene_symbol_report?hgnc_id=%s',
    EntrezGene  => 'https://www.ncbi.nlm.nih.gov/gene/%s', # Integer IDs
    EnsemblGene => 'http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s',
    OMIM       => 'https://omim.org/entry/%s',
    miRBase    => 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s',
    SnoRNAbase => 'https://www-snorna.biotoul.fr/plus.php?snoid=%s',
    MEROPS     => 'https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=%s',

    IMGT       => 'https://www.google.com/search?q=site:http://www.imgt.org/IMGTrepertoire/+%s',
    UCSC       => 'http://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=%s',
    MGD        => 'http://www.informatics.jax.org/marker/%s',
    RGD        => 'https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=%s',
};

my $typLU = {
    'withdrawn'           => 'Deprecated',
    'phenotype'           => 'Phenotype',
    'pseudogene'          => 'Pseudogene',
    'other'               => 'Other',
    'non-coding RNA'      => 'ncRNA',
    'protein-coding gene' => 'ProteinCoding',
};
my $taxa = {
    MGD => 'Mus musculus',
    RGD => 'Rattus norvegicus',
};

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
    ftp        => $defFtp,
    symlevels  => $defSym,
    typelevels => $defTyp,
    dir        => ".",
};

BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}
use lib "$scriptDir";
require Utils;
our ($args, $clobber, $ftp, $tmpDir, $maxAbst, $bar);

my $outDir   = $args->{dir};    $outDir =~ s/\/+$//;
my $stash    = $args->{stash};

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

      -ftp Default '$defFtp'. Domain of EBI's FTP site

  -clobber Default 0, which will preserve any already-generated
           files. Set to 1 to regenerate matrix files, and any value
           greater than 1 to also re-download base files from Entrez.

-symlevels Default '$defSym'.
           Factor levels for symbol nomenclature status.

-typelevels Default '$defTyp'.
           Factor levels for gene types.

    -stash If a non-zero argument, try to register the file in the Stash
           database, if available.

";
    exit;
}

my $src    = "pub/databases/genenames/new/tsv/hgnc_complete_set.txt";
my $srcUrl = join('/', "ftp:/", $args->{ftp}, $src);
my $loc    = &fetch_url($src);
my $vers   = &_datestamp_for_file($loc);

## File format is:
## <Type>@<RowNamespace>_to_<ColNamespace>@<Authority>@<Version>.mtx
##    <Type> = Type of Matrix, generally Map or Ontology
##    <RowNamespace> = namespace of the row identifiers
##    <ColNamespace> = namespace of the column identifiers
##    <Authority> = Entity providing the raw data
##    <Version> = Version token (eg GRCh37) or YYYY-MM-DD of source file(s)
my $fFmt   = sprintf('%s/%%s@%%s_to_%%s@HGNC@%s.mtx', $outDir, $vers);

warn "

Building pivot matrices from HGNC:

   FTP Site: $defFtp
     Source: $src
 Local Copy: $loc
  Datestamp: $vers

";

&symbol_maps();
&accession_maps();

sub accession_maps {
    # Grab the following namespaces:
    my ($rows, $xcol) =
        &_get_rows([ HGNC        => 'hgnc_id', # Easier to grab this value twice
                     OMIM        => 'omim_id',
                     miRBase     => 'mirbase',
                     MEROPS      => 'merops',
                     IMGT        => 'imgt',
                     UCSC        => 'ucsc_id',
                     MGD         => 'mgd_id',
                     RGD         => 'rgd_id',
                     EntrezGene  => 'entrez_id',
                     EnsemblGene => 'ensembl_gene_id' ]);

    ## Factorize the gene types
    my ($lvls, $scH) = &_factorize_levels('typelevels');

    my @gMetaCols = qw(Symbol Type Description);
    my $nsNum = $#{$xcol};
    ## Make all pairwise matrices
    for my $i (0..($nsNum-1)) {
        my $nsi = $xcol->[$i];
        for my $j (($i+1)..$nsNum) {
            my $nsj = $xcol->[$j];
            my $trg = sprintf($fFmt, "Map", $nsi, $nsj);
            my %tax = map { ($taxa->{$_} || "Homo sapiens") => 1 } ($nsi, $nsj);
            my $meta = {
                Species   => [sort keys %tax],
                Namespace => [$nsi, $nsj],
            };
            unless (&output_needs_creation($trg)) {
                ## File is already made, and clobber is not set
                &msg("Keeping existing accession map:", $trg);
                &stash($trg, $meta);
                next;
            }

            my (%data, %iMeta, %jMeta, %seenLvl);
            my ($jnum, $nznum) = (0, 0);
            foreach my $row (@{$rows}) {
                # Normalize the text by removing any quotes
                my ($iTxt, $jTxt) = ($row->[1][$i] || "", $row->[1][$j] || "");
                if ($iTxt =~ /^"\s*(.+?)\s*"$/ || $iTxt =~ /^'\s*(.+?)\s*'$/) {
                    $iTxt = $1;
                }
                if ($jTxt =~ /^"\s*(.+?)\s*"$/ || $jTxt =~ /^'\s*(.+?)\s*'$/) {
                    $jTxt = $1;
                }
                # Some genes are not represented in some namespaces:
                next unless ($iTxt && $jTxt);
                # Score will be based on the gene type
                my $type = $typLU->{ $row->[0][2] || "" } || "Unknown";
                my $sc   = $scH->{ uc($type) };
                $seenLvl{ $type }++;

                # Some maps will not be 1:1 for all IDs
                my @iids = split(/\s*\|\s*/, $iTxt);
                my @jids = split(/\s*\|\s*/, $jTxt);
                
                foreach my $iid (@iids) {
                    $iMeta{$iid} ||= {
                        Symbol => $row->[0][1] || $row->[0][4] || "",
                        Type   => $type,
                        Description => $row->[0][3],
                    };
                    my $targ =  $data{$iid} ||= {
                        name => $iid,
                        hits => [],
                    };
                    foreach my $jid (@jids) {
                        my $jdat = $jMeta{$jid} ||= {
                            Symbol => $row->[0][1] || $row->[0][4] || "",
                            Type   => $type,
                            Description => $row->[0][3],
                            num         => ++$jnum,
                        };
                        push @{$targ->{hits}}, ($jdat->{num}, $sc);
                        $nznum++;
                    }
                }
            }

            my @iids = sort keys %data;
            my @jids = sort { $jMeta{$a}{num} <=> $jMeta{$b}{num} } keys %jMeta;
            my ($rnum,$cnum) = ($#iids + 1, $#jids + 1);

            my $tmp = "$trg.tmp";
            open(MTX, ">$tmp") || &death("Failed to write Symbol map", $tmp, $!);
            print MTX &_initial_mtx_block
                ( "Mapping", $rnum, $cnum, $nznum, "HGNC $nsi-to-$nsj Map",
                  "Scores are factor levels representing the gene type.",
                  "$nsi ID", $nsUrl->{$nsi},
                  "$nsj ID", $nsUrl->{$nsj}, undef, $srcUrl);

            print MTX &_citation_MTX();
            print MTX &_species_MTX( $meta->{Species} );
            print MTX &_filter_block({
                KeepLevel => "[,][ncRNA,ProteinCoding] ## Most applications will want to work with coding and 'typical' ncRNA genes"
                                     });

            print MTX &_factor_map_block( $lvls, $scH, "Gene Type", undef,
                                          \%seenLvl);
            print MTX &_rowcol_meta_comment_block();

            ## Note row metadata
            print MTX &_generic_meta_block(\@iids, 'Row', \%iMeta, \@gMetaCols);

            ## Column metadata
            print MTX "% $bar\n";
            print MTX &_generic_meta_block(\@jids, 'Col', \%jMeta, \@gMetaCols);

            print MTX &_triple_header_block( $rnum, $cnum, $nznum );
            for my $i (0..$#iids) {
                my $hits = $data{$iids[$i]}{hits};
                for (my $j = 0; $j < $#{$hits}; $j += 2) {
                    printf(MTX "%d %d %d\n", $i+1, $hits->[$j], $hits->[$j+1]);
                }
            }
            close MTX;
            rename($tmp, $trg);
            &msg("Generated $nsi to $nsj mapping", $trg);
            &stash($trg, $meta);
            
        }
    }

}

sub symbol_maps {
    # Grab the following namespaces:
    my ($rows, $xcol) =
        &_get_rows([ HGNC        => 'hgnc_id', # Easier to grab this value twice
                     EntrezGene  => 'entrez_id',
                     EnsemblGene => 'ensembl_gene_id' ]);
    
    ## Factorize the symbols
    my ($lvls, $scH) = &_factorize_levels('symlevels');
    my $offL = $scH->{ uc("Official") };
    my $upL  = $scH->{ uc("UnofficialPreferred") };
    my $uL   = $scH->{ uc("Unofficial") };
    my $oldL = $scH->{ uc("Unknown") };
    foreach my $row (@{$rows}) {
        my %syms;
        ## Is the entry withdrawn? If so we will set all scores to Unknown:
        my $isWD = $row->[0][2] =~ /withdrawn/i ? $oldL : 0;
        if (my $off = $row->[0][1]) {$syms{$off} ||= $isWD || $offL;} # Official
        if (my $pre = $row->[0][4]) {$syms{$pre} ||= $isWD || $upL;}  # UnPref
        map { $syms{$_} ||= $isWD || $uL } @{$row->[0][5]}; # All other Unoff
        map { $syms{$_} ||= $oldL } @{$row->[0][6]}; # "Previous" -> Unknown
        $row->[2] = \%syms;
    }

    ## Now build the lookups for the 'main' namespaces
    my @gMetaCols = qw(Symbol Type Description);
    for my $xi (0..$#{$xcol}) {
        my $ns  = $xcol->[$xi];
        my $trg = sprintf($fFmt, "Map", "Symbol", $ns);
        my $meta = {
            Species   => ["Homo sapiens"],
            Namespace => ["Symbol", $ns],
        };
        unless (&output_needs_creation($trg)) {
            ## File is already made, and clobber is not set
            &msg("Keeping existing Symbol map:", $trg);
            &stash($trg, $meta);
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
              "Gene Symbol", $nsUrl->{Symbol},
              "$ns ID", $nsUrl->{$ns}, undef, $srcUrl);

        print MTX &_citation_MTX();
        print MTX &_species_MTX( ["Homo sapiens"] );
        print MTX &_filter_block({
            TossLevel => "[,][Unknown] ## 'Unknown' entries are deprecated assignments."
                                 });

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
                printf(MTX "%d %d %d\n", $i+1, $hits->[$j], $hits->[$j+1]);
            }
        }
        close MTX;
        rename($tmp, $trg);
        &msg("Generated Symbol to $ns mapping", $trg);
        &stash($trg, $meta);

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

sub _citation_MTX {
    return "% 
%% DEFAULT Citation HGNC Database, HUGO Gene Nomenclature Committee (HGNC), EMBL Outstation - Hinxton, European Bioinformatics Institute, Wellcome Trust Genome Campus, Hinxton, Cambridgeshire, CB10 1SD, UK www.genenames.org. Downloaded on $vers
%
";
}

sub _species_MTX {
    my ($tax) = @_;
    return $tax ?
        sprintf("%% DEFAULT Species [,][%s]\n",join(',',@{$tax})) : "";
}

sub stash {
    return unless ($stash);
    my $exe = `which stash`;
    unless ($exe) {
        warn "
stash is not installed on your system.
    Unless you were really expecting it to be there, this is not an error.

";
        $stash = 0; # Only warn once
        return;
    }
    my ($file, $mh) = @_;
    $mh ||= {};
    $mh->{MatrixType} ||= "Map",
    $mh->{Authority}    = "HGNC";
    $mh->{Version}      = $vers;
    $mh->{Format}       = "MatrixMarket";
    $mh->{FileType}     = "AnnotatedMatrix";
    my @meta;
    while (my ($k, $v) = each %{$mh}) {
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
    my $cmd = "$exe add --metadata '".join(',', @meta)."' \"$file\"";
    warn "$cmd\n";
}
