#!/usr/bin/perl -w

die "

I was unable to find current documentation for the BioMart restful
API. Abandoning this tool in favor of an R implementation using
'biomaRt'.

";

use strict;
my $scriptDir;
our $defTmp   = "/tmp/bioMartFiles";
my $defMart   = "http://www.ensembl.org/biomart/martservice";

our $defaultArgs = {
    dataset  => "",
    dir      => ".",
};

BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}
use lib "$scriptDir";
require Utils;
our ($args, $clobber, $ftp, $tmpDir, $maxAbst);

my ($datasets, $dsLU) = &gene_datasets();

my $dsDat = $dsLU->{lc($args->{dataset} || "")};


if ($args->{h} || $args->{help} || !($dsDat) ) {
    warn "
Usage:

This program will parse Ensembl information to generate lookup
matrices, written in MatrixMarket format.

Required Arguments:

  -dataset The Ensembl dataset to extract data from. Can be provided
           as the formal dataset name, the Ensembl common name, or the
           genome build name. A list of recognized datasets will be
           shown if a known one is not provided.

Optional Arguments:

   -tmpdir Default '$defTmp'. Directory holding downloaded files

      -dir Default '.'.  Output directory that will contain generated
           matrix and metadata files

  -clobber Default 0, which will preserve any already-generated
           files. Set to 1 to regenerate matrix files, and any value
           greater than 1 to also re-download base files from Entrez.

     -name A name to use for the files. By default this will be built
           from the 'common' name used by Ensembl for the dataset plus
           the genome build.

     -help Show this documentation

";
    unless ($dsDat) {
        my @bits = ("Please provide a -dataset. Available ones",
                    "  (can use any of the three values to specify):");
        push @bits, map { sprintf("%25s = %s (%s)", @{$_}) } (['<ORGANISM>','<DATASET>', '<BUILD>'], @{$datasets});
        &err(@bits);
    }
    exit;
}

my $outDir   = $args->{dir};    $outDir =~ s/\/+$//;
my $dataset  = $dsDat->[3];
my $filetok  = $args->{name};
unless ($filetok) {
    $filetok = join('-', $dsDat->[0], 'Ensembl', $dsDat->[2]);
    $filetok =~ s/\s+/_/g;
}

&mkpath([$outDir]);

&geneinfo_file();

sub geneinfo_file {
    my $trg = sprintf("%s/Metadata-%s_GeneInfo.tsv", $outDir, $filetok);
    unless (&output_needs_creation($trg)) {
        &msg("Using existing GeneInfo file:", $trg);
        return $trg;
    }
    &msg("Recovering GeneInfo file from BioMart");
    my $tmp = "$trg.tmp";
    my $qry = &biomart_geneinfo_query( $dataset );
    &get_url( $defMart, $trg, { query => $qry });
    &msg("GeneInfo file recovered", $trg);
    return $trg;
}

sub biomart_geneinfo_query {
    my $dataset = shift; # hsapiens_gene_ensembl
    my $rv = sprintf('
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1"
       uniqueRows="1" count="" datasetConfigVersion="0.6" >
    <Dataset name="%s" interface="default" >
        <Attribute name="ensembl_gene_id" />
        <Attribute name="ensembl_transcript_id" />
        <Attribute name="ensembl_peptide_id" />
        <Attribute name="description" />
        <Attribute name="external_gene_name" />
        <Attribute name="external_transcript_name" />
        <Attribute name="gene_biotype" />
        <Attribute name="transcript_biotype" />
        <Attribute name="status" />
        <Attribute name="transcript_status" />
    </Dataset>
</Query>
', $dataset);
    ## $rv =~ s/\s*[\n\r]+\s*//; # Remove newlines
    return $rv;
}

sub gene_datasets {
    # Have not yet found a query to allow these to be recovered
    # progamatically
    my $sourceScrape = '<option value="ggallus_gene_ensembl">Chicken genes (Gallus_gallus-5.0)</option><option value="hsapiens_gene_ensembl">Human genes (GRCh38.p10)</option><option value="mmusculus_gene_ensembl">Mouse genes (GRCm38.p5)</option><option value="rnorvegicus_gene_ensembl">Rat genes (Rnor_6.0)</option><option value="drerio_gene_ensembl">Zebrafish genes (GRCz10)</option><option value="">-------------------------------------------</option><option value="vpacos_gene_ensembl">Alpaca genes (vicPac1)</option><option value="pformosa_gene_ensembl">Amazon molly genes (Poecilia_formosa-5.1.2)</option><option value="acarolinensis_gene_ensembl">Anole lizard genes (AnoCar2.0)</option><option value="dnovemcinctus_gene_ensembl">Armadillo genes (Dasnov3.0)</option><option value="ogarnettii_gene_ensembl">Bushbaby genes (OtoGar3)</option><option value="cintestinalis_gene_ensembl">C.intestinalis genes (KH)</option><option value="csavignyi_gene_ensembl">C.savignyi genes (CSAV 2.0)</option><option value="celegans_gene_ensembl">Caenorhabditis elegans genes (WBcel235)</option><option value="fcatus_gene_ensembl">Cat genes (Felis_catus_6.2)</option><option value="amexicanus_gene_ensembl">Cave fish genes (AstMex102)</option><option value="ptroglodytes_gene_ensembl">Chimpanzee genes (CHIMP2.1.4)</option><option value="psinensis_gene_ensembl">Chinese softshell turtle genes (PelSin_1.0)</option><option value="gmorhua_gene_ensembl">Cod genes (gadMor1)</option><option value="lchalumnae_gene_ensembl">Coelacanth genes (LatCha1)</option><option value="btaurus_gene_ensembl">Cow genes (UMD3.1)</option><option value="cfamiliaris_gene_ensembl">Dog genes (CanFam3.1)</option><option value="ttruncatus_gene_ensembl">Dolphin genes (turTru1)</option><option value="aplatyrhynchos_gene_ensembl">Duck genes (BGI_duck_1.0)</option><option value="lafricana_gene_ensembl">Elephant genes (Loxafr3.0)</option><option value="mfuro_gene_ensembl">Ferret genes (MusPutFur1.0)</option><option value="falbicollis_gene_ensembl">Flycatcher genes (FicAlb_1.4)</option><option value="dmelanogaster_gene_ensembl">Fruitfly genes (BDGP6)</option><option value="trubripes_gene_ensembl">Fugu genes (FUGU 4.0)</option><option value="nleucogenys_gene_ensembl">Gibbon genes (Nleu1.0)</option><option value="ggorilla_gene_ensembl">Gorilla genes (gorGor3.1)</option><option value="cporcellus_gene_ensembl">Guinea Pig genes (cavPor3)</option><option value="eeuropaeus_gene_ensembl">Hedgehog genes (eriEur1)</option><option value="ecaballus_gene_ensembl">Horse genes (Equ Cab 2)</option><option value="pcapensis_gene_ensembl">Hyrax genes (proCap1)</option><option value="dordii_gene_ensembl">Kangaroo rat genes (dipOrd1)</option><option value="pmarinus_gene_ensembl">Lamprey genes (Pmarinus_7.0)</option><option value="etelfairi_gene_ensembl">Lesser hedgehog tenrec genes (TENREC)</option><option value="mmulatta_gene_ensembl">Macaque genes (Mmul_8.0.1)</option><option value="cjacchus_gene_ensembl">Marmoset genes (C_jacchus3.2.1)</option><option value="olatipes_gene_ensembl">Medaka genes (HdrR)</option><option value="pvampyrus_gene_ensembl">Megabat genes (pteVam1)</option><option value="mlucifugus_gene_ensembl">Microbat genes (Myoluc2.0)</option><option value="mmurinus_gene_ensembl">Mouse Lemur genes (Mmur_2.0)</option><option value="panubis_gene_ensembl">Olive baboon genes (PapAnu2.0)</option><option value="mdomestica_gene_ensembl">Opossum genes (monDom5)</option><option value="pabelii_gene_ensembl">Orangutan genes (PPYG2)</option><option value="amelanoleuca_gene_ensembl">Panda genes (ailMel1)</option><option value="sscrofa_gene_ensembl">Pig genes (Sscrofa10.2)</option><option value="oprinceps_gene_ensembl">Pika genes (OchPri2.0-Ens)</option><option value="xmaculatus_gene_ensembl">Platyfish genes (Xipmac4.4.2)</option><option value="oanatinus_gene_ensembl">Platypus genes (OANA5)</option><option value="ocuniculus_gene_ensembl">Rabbit genes (OryCun2.0)</option><option value="scerevisiae_gene_ensembl">Saccharomyces cerevisiae genes (R64-1-1)</option><option value="oaries_gene_ensembl">Sheep genes (Oar_v3.1)</option><option value="saraneus_gene_ensembl">Shrew genes (sorAra1)</option><option value="choffmanni_gene_ensembl">Sloth genes (choHof1)</option><option value="loculatus_gene_ensembl">Spotted gar genes (LepOcu1)</option><option value="itridecemlineatus_gene_ensembl">Squirrel genes (spetri2)</option><option value="gaculeatus_gene_ensembl">Stickleback genes (BROAD S1)</option><option value="csyrichta_gene_ensembl">Tarsier genes (tarSyr1)</option><option value="sharrisii_gene_ensembl">Tasmanian devil genes (Devil_ref v7.0)</option><option value="tnigroviridis_gene_ensembl">Tetraodon genes (TETRAODON 8.0)</option><option value="oniloticus_gene_ensembl">Tilapia genes (Orenil1.0)</option><option value="tbelangeri_gene_ensembl">Tree Shrew genes (tupBel1)</option><option value="mgallopavo_gene_ensembl">Turkey genes (Turkey_2.01)</option><option value="csabaeus_gene_ensembl">Vervet-AGM genes (ChlSab1.1)</option><option value="neugenii_gene_ensembl">Wallaby genes (Meug_1.0)</option><option value="xtropicalis_gene_ensembl">Xenopus genes (JGI 4.2)</option><option value="tguttata_gene_ensembl">Zebra Finch genes (taeGut3.2.4)</option>';
    my ($sets, $lu) = @_;
    foreach my $bit (split(/<\/option>/, $sourceScrape)) {
        if ($bit =~ /value="(.+?)">(.+) genes \((.+?)\)/) {
            my ($dataset, $name, $bld) = ($1, $2, $3);
            ## OMG build subversions have almost no standardization
            ## GRCh38        -> GRCh38.p10
            ## Gallus_gallus -> Gallus_gallus-5.0
            ## C_jacchus     -> C_jacchus3.2.1
            ## Stickleback   -> 'BROAD S1'
            ## Xenopus       -> 'JGI 4.2'
            my $dsmall = $dataset; $dsmall =~ s/_.+//;
            my $dat = [$name, $dsmall, $bld, $dataset];
            push @{$sets}, $dat;
            map { $lu->{lc($_)} = $dat } @{$dat};
            
        }
    }
    return ($sets, $lu);
}
