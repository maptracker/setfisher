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
our ($args, $clobber, $ftp, $tmpDir, $maxAbst, $bar);
my (%doneStuff);

my $revision = '1';
my $revNotes = {
    '0' => "# Beta code, still under development",
    '1' => "# Functional matrices being generated",
};

use XML::Parser::PerlSAX;
# print Dumper($args); die;

my ($sec, $min, $hr, $day, $mon, $year, $wday) = localtime;
my $today = sprintf("%04d-%02d-%02d", $year+1900, $mon+1, $day);

my $xml     = $args->{xml};


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
        &err("Please provide the path to the XML file");
        warn "
  XML files can be freely downloaded here:

    http://software.broadinstitute.org/gsea/downloads.jsp#msigdb

  You should download 'All gene sets' in 'Current MSigDB xml file' format
  
  The file name should look something like 'msigdb_v6.2.xml'. You will need 
  to register with the Broad before doing this
  (which is why the script can't perform the download for you).

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
    PMID             => 'SourcePMID',
    AUTHORS          => 'SourceAuthors',
    GEOID            => 'GEOID',
    EXACT_SOURCE     => 'SourceDetails',
    GENESET_LISTING_URL => 'SourceURL',
    CHIP                => 'World',
    SUB_CATEGORY_CODE   => 'Subcategory',
    DESCRIPTION_BRIEF   => 'Description',
    DESCRIPTION_FULL    => 'DescriptionFull',
    ORGANISM            => 'Species',
};

my @stndRowMeta = qw(Symbol OtherNames);
my @stndColMeta = qw(Subcategory Species World Description GEOID SourcePMID DescriptionFull SourceDetails SourceURL SourceAuthors OldNames);

my $defColDef = {
    ## Gene metadata
    Symbol => "Human Entrez Gene symbol, as normalized / chosen by MSigDB",
    OtherNames => "Other symbols or gene names found to be associated with the gene, including those from orthologues in other species",

    ## Gene set metadata
    Subcategory => 'MSigDB subcategory, generally an external database',
    Species => "The organism the set was found in. Note that the set member Entrez Gene IDs are all human, however",
    World => 'Usually a chip identifier, an indication of the gene space the list was taken from',
    Description => 'The brief MSigDB description for the gene set',
    GEOID => 'GEO or ArrayExpress ID',
    SourcePMID => 'The PubMed publication ID, if available',
    DescriptionFull => 'The full MSigDB description, may be a paper abstract',
    SourceAuthors => 'Authors on the source publication',
    SourceDetails => 'Additional information, such as a figure or table number, detailing how the list was taken from the source',
    SourceURL => 'The URL to the source',
    OldNames => 'Gene set names for this set from older versions of MSigDB',
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

## Hyperlink templates
my $nsUrl = {
    EntrezGene     => 'https://www.ncbi.nlm.nih.gov/gene/%s', # Integer IDs
    MSigDB         => 'http://software.broadinstitute.org/gsea/msigdb/cards/%s',
};

## Top-level hash to organize the data collected by the SAX parser:
my $saxData = { };

&parseXML();

## Static metadata
my $nsi       = "EntrezGene"; # Row namespace
my $nsj       = "MSigDB";     # Col namespace
my $type      = "Ontology";   # Matrix type
my $auth      = "MSigDB";     # Primary data authority
my $authLong  = "$auth ## Initiative at the Broad Institute to provide currated gene sets and tools to analyze them";

## Metadata dependent on the particular XML file being parsed
my $vers      = $saxData->{build}{Version};
my $versToken = $vers;
my $source    = sprintf("http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/%s/msigdb_v%s.xml",
 $versToken, $versToken);

## Stash deduplicated file store - not on all systems
my $stashMeta = {
    Authority  => $auth,
    'Version'  => $vers,
    MatrixType => $type,
    Revision   => $revision,
    FileType   => "AnnotatedMatrix",
    'Format'   => "MatrixMarket",
};



&generateMatrices();

&msg("  Done: All matrices created");

sub generateMatrices {
    foreach my $cH (sort {$a->{Category} cmp $b->{Category} }
                    values %{$saxData->{sets}}) {
        &generateMatrix( $cH );
    }
}

sub generateMatrix {
    my $cH      = shift;
    my $cat     = $cH->{Category};
    my $catName = $cH->{CategoryName};
    my $mod     = sprintf("%s.%s", $cat, $catName);
    my %fbits   = (type => $type,  mod  => $mod,
                   ns1  => $nsi,   ns2  => $nsj,
                   auth => $auth,  vers => $versToken,  
                   dir => "$auth/$versToken");
    
    my $meta = {
        MatrixType => $fbits{type},
        Modifier   => $fbits{mod},
        Source     => $source,
        Namespace  => [$nsi, $nsj],
    };
    my $fmeta = { %{$stashMeta}, %{$meta} };
    my $trg   = &primary_path(%fbits);
    unless (&output_needs_creation($trg)) {
        &msg("Keeping existing $nsj file:", $trg);
        &post_process( %fbits, meta => $fmeta );
        return;
    }

    my %rm     = %{$cH->{genes}}; # Row meta
    my %cm     = %{$cH->{sets}};  # Col meta
    my @rids   = sort { $rm{$a}{order} <=> $rm{$b}{order} } keys %rm;
    my @cids   = sort { $cm{$a}{order} <=> $cm{$b}{order} } keys %cm;

    foreach my $rH (values %rm) {
        ## Normalize the OtherNames metadata for genes
        my %on = %{$rH->{other}};
        delete $on{ $rH->{Symbol} || ""}; # Remove the 'main' symbol from names
        $rH->{OtherNames} = join(',', sort {uc($a) cmp uc($b) } keys %on) || "";
    }
    
    my ($rnum, $cnum, $nznum) = ($#rids + 1, $#cids + 1, $cH->{nznum});
    
    my $tmp = "$trg.tmp";
    open(FH, ">$tmp") || &death("Failed to write $nsi $nsj $type", $tmp, $!);
    print FH &_initial_mtx_block
        ( $type, $rnum, $cnum, $nznum, "$nsi to $mod $nsj $type",
          "$auth collection $cat ($catName), $cH->{Description}",
          "Scores are all 1, no additional information on strength of assignment are available",
          $nsi, $nsj);

    print FH &_dim_block({
        %{$fmeta},
        RowDim    => $nsi,
        RowUrl    => $nsUrl->{$nsi},
        ColDim    => $nsj,
        ColUrl    => $nsUrl->{$auth}, 
        Authority => $authLong },
        { Version  => "# $saxData->{build}{BuildDate}",
          Revision => $revNotes->{defined $revision ? $revision : ""} } );
    
    ## Some basic information on the set
    foreach my $setInfo (qw(Category CategoryName)) {
        print FH &_default_parameter($setInfo, $cH->{$setInfo});
    }

    print FH "%\n";
    print FH &_citation_MTX();

    ## Just include populated metadata for the columns. Some
    ## collections (eg C1.Positional) have sparser annotation than
    ## others.
    my @usedColMeta;
    foreach my $scm (@stndColMeta) {
        for my $c (0..$#cids) {
            if ($cm{$cids[$c]}{$scm}) {
                push @usedColMeta, $scm;
                last;
            }
        }
    }

    my $usedColDef = { map {$_ => $defColDef->{$_} } 
                       (@stndRowMeta, @usedColMeta) };

    print FH &_rowcol_meta_comment_block( $usedColDef );
    ## Row metadata
    print FH &_generic_meta_block(\@rids, 'Row', \%rm, \@stndRowMeta);
    print FH "% $bar\n";

    print FH &_generic_meta_block(\@cids, 'Col', \%cm, \@usedColMeta);
    print FH "% $bar\n";
    
    ## Make the Row-Col connections. All scores are '1'
    print FH &_triple_header_block( $rnum, $cnum, $nznum );
    foreach my $rid (@rids) {
        my $gdat = $rm{$rid};
        my $rInd = $gdat->{order};
        foreach my $cInd (sort {$a <=> $b} keys %{$gdat->{sets}}) {
            printf(FH "%d %d 1\n", $rInd, $cInd);
        }
    }
    close FH;
    rename($tmp, $trg);
    &post_process( %fbits, meta => $fmeta );
    &msg("Generated $mod $nsi -> $nsj", $trg);
}


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


sub _citation_MTX {
    return &_default_parameter( "Citation", "Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles; Subramanian et al, PNAS Oct 25, 2005. 102 (43) 15545-15550 / Molecular signatures database (MSigDB) 3.0; Liberzon et al, Bioinformatics, Volume 27, Issue 12, 15 Jun 2011, Pages 1739–1740 / The Molecular Signatures Database Hallmark Gene Set Collection; Liberzon et al, Cell Systems, Volume 1, Issue 6, 23 Dec 2015, Pages 417-425 / $source")."%\n";
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
    
    ## Dangit dangit dangit. I believed that the species indicated the
    ## taxa of the Entrez GeneID. It appears that all gene IDs are,
    ## however, human. This is generally good, in the sense that it
    ## normalizes gene IDs and allows set information to be pulled
    ## from biological data in other species. It means we shouldn't
    ## segregate the gene sets by species, though, and instead the
    ## species should be a property of each set.

    ## my $catSpec = "$cat\t$spec"; # nope
    my $catSpec = $cat; # Just the category
    
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
            Description  => $cD,
            geneCount    => 0,
            setCount     => 0,
            sets         => {},
            genes        => {},
            nznum        => 0,
        };
    }
    my $cH    = $data->{sets}{$catSpec};
    my $set   = $attr{STANDARD_NAME} || 'ERROR_NO_SET_NAME';
    if ($cH->{sets}{$set}) {
        warn "$cat set $set is apparently present twice\n";
        return;
    }
    my $sdat = $cH->{sets}{$set} = {
        order => ++$cH->{setCount}, # Matrix index for this set
        noMap => {}, # Member names that could not be mapped to Entrez
    };
    my $sNum = $sdat->{order};
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
                order  => ++$cH->{geneCount}, ## Matrix index for this gene
                Symbol => $sym,
                other  => {},
                sets   => {},
            };
            $gdat->{other}{$n}++ if ($n);
            unless ($gdat->{sets}{$sNum}++) {
                ## Tally total number of unique connections
                $cH->{nznum}++;
            }
        } else {
            $sdat->{noMap}{$n}++;
        }
    }
}
