#!/usr/bin/perl -w

use strict;
my $scriptDir;
our $defTmp   = "/tmp/msigdbFiles";

our $defaultArgs = {
};

BEGIN {
    use File::Basename;
    $scriptDir = dirname($0);
}
use lib "$scriptDir";
require Utils;
our ($args, $clobber, $ftp, $tmpDir, $maxAbst);
my (%doneStuff);

use XML::Parser::PerlSAX;
# print Dumper($args); die;

my ($sec, $min, $hr, $day, $mon, $year, $wday) = localtime;
my $today = sprintf("%04d-%02d-%02d", $year+1900, $mon+1, $day);


my $xml     = $args->{xml};
my $fileReq = $args->{file};
my $pmdb    = $args->{pubmeddb};
my $email   = $args->{email};
my $analyze = $args->{analyze};
my $keepxml = $args->{keepxml};

if ($args->{h} || $args->{help} || !$xml ) {
    warn "
Usage:

This program will parse a Broad MSigDB XML file to generate a set of
AnnotatedMatrix files.

Required Arguments:

      -dir Output directory that will contain generated matrix and
           metadata files

      -xml Path to MSigDB XML file. These files are freely available,
           but require registration to access.

Optional Arguments:

    -stash If a non-zero argument, try to register the file in the Stash
           database, if available.

     -help Show this documentation

";
    unless ($xml) {
        &err("Please provide the path to the XML file") unless ($xml);
        warn "
  XML files can be freely downloaded here:

    http://software.broadinstitute.org/gsea/downloads.jsp

  You should download 'All gene sets' in 'Current MSigDB xml file'
  format. The file name should look something like
  'msigdb_v6.2.xml'. You will need to register with the Broad before
  doing this (which is why the script can't perform the download for
  you).

";
    }
    exit;
}

&death("Provided XML file path does not exist or is empty") unless (-s $xml);

## Collection names and descriptions:
##   http://software.broadinstitute.org/gsea/msigdb/
my $collectMeta = {
    'H' => ['Hallmark', 'Coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.'],
        'C1' => ['Positional', 'Sets for each human chromosome and cytogenetic band.'],
        'C2' => ['Curated', 'Online pathway databases, publications in PubMed, and knowledge of domain experts.'],
        'C3' => ['Motif', 'Based on conserved cis-regulatory motifs from a comparative analysis of the human, mouse, rat, and dog genomes.'],
        'C4' => ['Computational', 'Defined by mining large collections of cancer-oriented microarray data.'],
        'C5' => ['GeneOntology', 'Genes annotated by the same GO terms.'],
        'C6' => ['Oncogenic', 'Defined directly from microarray gene expression data from cancer gene perturbations.'],
        'C7' => ['Immunologic', 'Defined directly from microarray gene expression data from immunologic studies.'],
};


## Attribute information:
##    https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_XML_description
my $attMap = {
    HISTORICAL_NAMES => 'OldNames',
    ORGANISM         => 'Species',
    PMID             => 'SourcePMID',
    AUTHORS          => 'SourceAuthors',
    GEOID            => 'GEOID',
    EXACT_SOURCE     => 'SourceDetails',
    GENESET_LISTING_URL => 'SourceURL',
    CHIP                => 'World',
    SUB_CATEGORY_CODE   => 'Subcategory',
    DESCRIPTION_BRIEF   => 'Description',
    DESCRIPTION_FULL    => 'DescriptionFull',
};

my $attDesc = {
    'OldNames' => 'Gene set names for this set from older versions of MSigDB',
'Species' => 'The organism the list was identified in',
'SourcePMID' => 'The PubMed publication ID, if available',
'SourceAuthors' => 'Authors on the source publication',
'GEOID' => 'GEO or ArrayExpress ID',
'SourceDetails' => 'Additional information, such as a figure or table number, detailing how the list was taken from the source',
'SourceURL' => 'The URL to the source',
'World' => 'Usually a chip identifier, an indication of the gene space the list was taken from',
'Subcategory' => 'MSigDB subcategory',
'Description' => 'The brief MSigDB description for the gene set',
'DescriptionFull' => 'The full MSigDB description, may be a paper abstract',
};


my $commonNames = {        # Counts from vers 6.2 - Jul 12, 2018
    "Homo sapiens"      => 'Human',         # 14660
    "Mus musculus"      => 'House mouse',   # 3974
    'Rattus norvegicus' => 'Norway rat',    # 29
    'Macaca mulatta'    => 'Rhesus monkey', # 6
    'Danio rerio'       => 'Zebrafish'      # 5
};


## Expand the subcategories to be more immediately interpretable
##   http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp
my $subCatNames = {
        'CC' => "Cellular Compartment",
        'CP:BIOCARTA' => "BioCarta",
        'CM' => "Cancer Modules",
        'MF' => "Molecular Function",
        'CP:REACTOME' => "Reactome",
        'CP:KEGG' => "KEGG",
        'CP' => "Canonical Pathways",
        'CGN' => "Cancer Gene Neighborhoods",
        'CGP' => "Chemical and Genetic Pertubations",
        'MIR' => "miRBase",
        'TFT' => "Transcription Factor Targets",
        'BP' => "Biological Process"
};

## Top-level hash to organize the data collected by the SAX parser:
my $saxData = { };

&parseXML();
print Dumper($saxData->{subCat}); die;


sub parseXML {
    my $t = time;
    my $handler = MSigDBHandler->new( $saxData );
    my $parser  = XML::Parser::PerlSAX->new( Handler => $handler );

    #eval {
        $parser->parse( Source => {SystemId => $xml } );
    #};
    if ($handler->{doneXML}) {
        my @bits = ("    XML: $xml",
                    "Records: $handler->{count}",
                    sprintf("   Time: %.1f min", (time - $t) / 60));
        &msg(@bits);
    } else {
        &death("Failed to complete parse of XML", $xml);
    }
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###   Custom SAX handler module
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

=head1 MSigDB XML Handler

Parsing is simplified in that all information are stored as XML
attributes.

=cut

package MSigDBHandler;

use Data::Dumper;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        count   => 0,
    };
    bless ($self, $class);
    my $data = $self->{data} = shift;
    return $self;
}

sub start_element {
    my ($self, $element) = @_;
    if ($element->{Name} eq 'GENESET') {
        $self->process_set( $element );
        $self->{count}++;
    } elsif ($element->{Name} eq 'MSIGDB') {
        ## Just a couple tidbits of information for the MSigDB release:
        $self->{data}{build}{Version}   = $element->{Attributes}{VERSION};
        $self->{data}{build}{BuildDate} = $element->{Attributes}{BUILD_DATE};
    }
}

sub end_element {
    my ($self, $element) = @_;
    if ($element->{Name} eq 'MSIGDB') {
        ## Just note that we made it to the end of the document
        $self->{doneXML} = 1;
    }
}

sub process_set {
    my $self = shift;
    my $node = shift;
    my $data = $self->{data};
    my %attr = %{$node->{Attributes}};
    my $cat  = $attr{CATEGORY_CODE} || 'UNK';
    my $spec = $attr{ORGANISM} || "SPECIES NOT SET";
    if ($cat eq 'ARCHIVED') {
        ## Older gene sets, looks like all are from deprecated GO terms
        my $subCat = $attr{SUB_CATEGORY_CODE} || "UNK";
        $data->{archived}{$subCat}++;
        return;
    }
    
    ## Tally up subcategories to see what they might represent:
    $data->{subCat}{$attr{SUB_CATEGORY_CODE} || "UNK"}++;
    ## Similarly look at species 
    $data->{specCount}{$spec}++;
    
    my $gbCom = $commonNames->{$spec};
    unless ($gbCom) {
        ## There are only a handful of species represented. Maintain
        ## an explicit lookup hash to recover common name.
        warn "Unrecognized species '$spec'
  Please update \$commonNames with the Genbank Common Name for this organism\n"
            unless ($doneStuff{"Species $spec"}++);
        return;
    }
    my $catSpec = "$cat\t$spec";
    unless ($data->{sets}{$catSpec}) {
        ## Data structure to aggregate information for each
        ## matrix. There will be one matrix per Category (aka
        ## Collection) / species.
        my ($cN, $cD) = @{$collectMeta->{$cat} || []};
        unless ($cN) {
            ## Insist that we know what each collection is:
            warn "No information available for Collection '$cat'
  Please update \$collectMeta to record information for this category.\n"
                unless ($doneStuff{"Category $cat"}++);
            return;
        }
        
        $data->{sets}{$catSpec} ||= {
            Category     => $cat,
            CategoryName => $cN,
            Species      => $spec,
            CommonName   => $gbCom,
            Description  => $cD,
            geneCount    => 0,
            setCount     => 0,
            sets         => {},
            genes        => {},
        };
    }
    my $cH    = $data->{sets}{$catSpec};
    my $set   = $attr{STANDARD_NAME} || 'ERROR_NO_SET_NAME';
    if ($cH->{sets}{$set}) {
        warn "$cat set $set in $gbCom is apparently present twice\n";
        return;
    }
    my $sdat = $cH->{sets}{$set} = {
        num   => ++$self->{setCount}, ## Matrix index for this set
        noMap => {}, ## Member names that could not be mapped to Entrez
    };
    ## Basic metadata for the set:
    while (my ($in, $out) = each %{$attMap}) {
        if (my $val = $attr{$in}) {
            ## We are just mapping over the attribute name here:
            $sdat->{$out} = $val;
        }
    }
    if (my $sc = $sdat->{Subcategory}) {
        if (my $nice = $subCatNames->{$sc}) {
            $sdat->{Subcategory} = $nice;
        } else {
            warn "Unrecognized subcategory '$sc'
  Please update \$subCatNames"
                unless ($doneStuff{"subcat $sc"}++);
        }
    }

    foreach my $mems (split(/\|/, $attr{MEMBERS_MAPPING} || "")) {
        my ($n, $sym, $gid) = split(',', $mems);
        if ($gid) {
            my $gdat = $cH->{genes}{$gid} ||= {
                num    => ++$self->{geneCount}, ## Matrix index for this gene
                Symbol => $sym,
                other  => {},
            };
            $gdat->{other}{$n}++;
        } else {
            $sdat->{noMap}{$n}++;
        }
    }
}
