### Background

Like most enrichment tools, SetFisher will ultimately be calculating
statistics using Fisher's Exact Test (the hypergeometric
distribution), via the `phyper()` method in R. This is a
"marble-counting" exercise, where a handful of black or white marbles
are pulled from an urn with a fixed, finite number, and the relative
distribtion of black and white in our handful is compared to the
distribution of the full urn. This generalized system can be applied
to many applications, but for us the situation will be:

* The marbles are genes
* The urn is the genome, or a restricted set of genes defined by a
  focused technology or protocol that is only assaying some of the
  genes (eg an Affymetrix array design)
* The handful of selected marbles represent "interesting" genes
  identified through an experiment; In many cases these will be genes
  that have undergone significant changes between a control sample and
  a treatment.
* The colors are assignments of the genes to a particular ontology,
  such as "localized to cytoplasm" or "cell death inhibitor"

In any Fisher's test, four integer values need to be calculated for
each test. In our case, a test is a pair between one gene list and one
ontology term, and the four values are:

1. `W` (sometimes refered to as `n+m`), the total number of genes in
   the "world" (the genome, the array design, etc).
1. `n`, the number of genes in the world assigned to the ontology term
1. `N`, the number of genes in selected list
1. `i`, the number of genes in the selection assigned to the ontology term

In R's phyper implementation, the arguments above will map as (R =
SetFisher):

1. "Selected white balls":
    * `q` = `i` + `lower.tail = TRUE` (under-enriched) _OR_
    * `q` = `i-1` + `lower.tail = FALSE` (overenriched)
1. "White balls in urn": `m` = `n`
1. "Black balls in urn": `n` = `W - n`
1. "Size of selection": `k` = `N`

### Terminology

_The terms below are explained in the context of genes, but the system
is agnostic and can be applied to any collection of objects and
ontologies._

* __Accession__ - a unique, reliable identifier that can be used and
  shared as an unambiguous ID for an "object" (gene, ontology term,
  publication, reagent, etc).
* __Description__ - human-readable "free text" that describes an
  object. Typically associated with __Accessions__.
* __Namespace__ - a source providing __Accessions__. Examples include:
    * Entrez Gene, eg "LOC859"
    * RefSeq Protein, eg "NP_001225"
    * Ensembl Gene, eg "ENSG00000182533"
    * Ensembl RNA, eg "ENST00000397368"
    * Gene Ontology, eg "GO:0072584"
    * Affymetrix Probeset, eg "41072_at"
* __Ontology Term__ - an __Accession__ associated with a specific
  property, behavior, classification, etc. The Ontology Term will have
  zero or more __IDs__ (genes in our case) assigned to it. Examples of
  gene-centric Ontology terms include:
    * "meiotic mismatch repair"
    * "WHEP-TRS domain"
    * "Crohn's Disease and Ulcerative Colitis"
* __Ontology__ - a _collection_ of __Ontology Terms__. Examples include:
    * Gene Ontology Biological Process
    * MSigDB Chromosomal Location
    * Wikipathways
    * PubMed (each article can have zero or more genes assigned to it)
* __ID__ - An __Accession__ utilized in the analysis. The ID
  __Namespace__ is always used for building __Ontology__
  matrices. __Queries__ may or may not be in the ID Namespace; If not,
  then a __Mapping Matrix__ is used to connect the __Query ID__
  Namespace to the ID Namespace.
* __Query ID__ - A gene identifier as provided in the input
  __Queries__. These identifiers _might_ be the same __Namespace__ as
  the IDs; If not, then a __Mapping Matrix__ will be needed to map
  Query IDs to IDs.
* __Mapping Matrix__ - a quantified lookup table mapping each __Query
  ID__ to zero or more __IDs__. See below for more details.
* __Query__ - sometimes refered to as "Lists" or "Gene Lists", a set
  of "interesting" __Query IDs__ derived from some observation or
  experiment. Queries have __List Names__ (either user-provided or
  auto-generated) assigned to them; They are _not_ referenced by
  __Accession__.
* __Analysis__ - A combination of one or more __Queries__ with one or
  more __Ontologies__, as well as a __Mapping Matrix__ if the __Query
  IDs__ have a different __Namespace__ than the gene __IDs__ used in
  the __Ontology__.

### Inputs

SetFisher is utilizing sparse matrices from the `Matrix` package to
encode set data. For any set enrichment implementation two classes of
sets (the queries and the ontologies) _must_ be provided, while a
third structure (the __Mapping Matrix__) may be requried in some
cases. As matrices these appear as:

1. One or more sets of interesting genes, the __Queries__, supplied as
   / stored in the __Query Matrix__:
    * Gene __Query IDs__ are in rows and __List Names__ are in columns
    * Non-zero cells indicate assignment of the row's gene to the
      column's list
    * Cell values could further be used for user-controlled
      thresholding _(Not yet implemented)_
1. One or more gene sets classifiers, the __Ontology__, stored as the
   __Ontology Matrix__:
    * Gene __IDs__ are in rows, __Ontology Terms__ in columns
    * Non-zero cells represent an assignment between the ID and
      Ontology Term.
    * The numeric value given to an assignment can be used to filter
      out some assignments
1. An optional __Mapping Matrix__ for __Namespace__ conversion
    * Maps the __Query ID__ namespace to the __ID__ namespace
    * Non-zero cells indicate that the Query ID can be mapped to the
      corresponding ID
    * Cell values can represent confidence scores and used for
      filtering out some mappings.

#### Why is a mapping matrix used?

The __Mapping Matrix__ could be ignored completely as long as
__Ontologies__ are built in the same __Namespace__ as the __Queries__
that are being used. The system could be used in such a way, but it is
not recommended for the following reasons:

* The mapping matrix allows abstraction of ontologies, which allows
  for significant efficiencies in large systems
    * For example, a single "Human Entrez GeneOntology" matrix can be
      created linking GeneOntology terms to NCBI Entrez Gene IDs. If
      an RNAseq experiment is reporting gene lists as Ensembl
      transcript IDs, rather than regenerating the ontology linked to
      ENST IDs, an "Ensemble RNA -> Entrez Gene" mapping matrix can be
      used
* Confidence scores captured in the matrix can be used to dynamically
  choose different assignment thresholds during analysis without
  requiring updates to the matrices.
    * For example, our Affymetrix-to-Gene mapping matrices record the
      highest fraction of Probeset probes that perfectly align to any
      of the gene's transcript. We normally consider a gene to hit if
      at least 880% of the probesets align to the gene (>= 0.8). That
      threshold is arbitrary, and can be relaxed or made more
      stringent if the user desires without having to regenerate
      either the mapping or ontology matrices.

But the primary reason for making the matrix is:

* The Mapping Matrix can be used to manage 1:many and many:many
  mappings between the __Query IDs__ and ontology __IDs__.
    * Not all __Namespaces__ are appropriate for Fisher's Exact Test!!
      The null hypothesis in phyper presumes independence of the
      "marbles". In biological sciences this is never perfectly true,
      but some situations are much less perfect than others. In
      particular, some Namespaces have rampant "multiple voting"
      problems:
        * RNA / Protein - some genes have only a single splice variant
          (one vote), while others have dozens.
        * RNA- or protein-based probes - (eg Affymetrix) some genes
          have a single probe associated with them, others have many.


### Pruning and Transforming Matrices

Once all matrices are available, additional steps are applied that
serve to reduce the size and complexity of the matrices, and in the
case of mapping matrices generate non-integer intermediate matrices to
address multiple voting issues. There are several user-configurable
parameters that affect this process:


* __Query Matrix__ filters are applied once when a query list is provided:
    * __minQueryScore__ : Minimum score requried to keep a __Query
      ID__ in a query list
    * __maxQueryScore__ : Maximum score allowed to keep a __Query ID__
      in a query list (for example if the list is quantified by a
      p-value)
* __Mapping Matrix__ fitlers are applied once to filter out
  low-scoring __Query ID__ to __ID__ mapping assignments.
    * __minMapMatch__ : Minimum score required to keep a map matrix
      assignment. Any scores that fall below this value will be
      removed, effectively "breaking" the link between the __Query
      ID__ and the __ID__.
* __Ontology Matrix__ filters are applied recursively via
  `.pruneOntologyMatrix()` and `.pruneMappingMatrix()`, and will
  repeat as long as prior rounds of filtering have resulted in removal
  of at least one __ID__ or __Ontology Term__. 
    * __maxSetPerc__ : Maximum % of world that can be assigned to a
      term. Any __Ontology Term__ that has been assigned more than
      this percentage of __IDs__ will be removed from the __Ontology__
      matrix. Designed to eliminate (often useless) generic ontology
      terms like 'Cell' or 'Chromosome'.
    * __minSetSize__ : Minimum count of set members required in a
      term. Any __Ontology Term__ that has fewer than this number (not
      percent) of __IDs__ will be removed from the
      __Ontology__. Removes poorly-populated terms (think "widget
      transporting protein, type IX, subtype B, lime flavored" with 1
      gene assigned)
    * __minOntoSize__ : Minimum number of terms a set member should
      have. Removes gene __IDs__ that have fewer than this number of
      __Ontology Terms__. Intended to remove speculative,
      poorly-annotated "marbles" from the urn. __This parameter will
      likely have the highest impact on most analyeses__. Statistics
      will often be skewed by the presence of "fanciful" genes that
      are either rarely expressed or simply fictional. These
      "semi-genes" are akin to marbles glued to the bottom of the
      urn - _you can never select them_. This has the effect of
      inflating `N`, which tends to inflate significance (lower
      p-values).
    * __minOntoMatch__ : Minimum score requried to keep an ontology
      matrix assignment. This is a simple assignment filter that sets
      cells in the __Ontology Matrix__ that fall below the threshold
      to zero (ie, breaks the link between the __ID__ and the
      __Ontology Term__).

#### Consequences of pruning

Effectively all parts of the contingency table are affected by pruning:

* `W` (surviving __ID__ world size) will vary by ontology! If an
  ontology provides "poor support" for some genes (few __Terms__
  assigned to the relevent __IDs__), then __minOntoSize__ can result
  in a large number of __IDs__ being pruned from the world
* The number of terms in an ontology will vary according to
  world. Smaller worlds will make it more likely that ontology terms
  will fail a __minSetSize__ filter.

### Fractional Counting

If a __Mapping Matrix__ is not utilized, then __Query IDs__ are the
same as __IDs__, and every __Query ID__ will "count as one marble" -
it will contribute a discrete count in the tallies of `W`, `N`, `n`
and `i`.

In other cases when a mapping matrix is used, discrete counting will
still be utilized. For example, the Query IDs are official gene
symbols, the ontology is represented with Entrez Gene IDs, and the
mapping matrix is simply a 1:1 "lookup" that translates each symbol to
the appropriate ID.

However, the primary impetus for implementing the mapping matrix was
to manage many:many relationships between the Query IDs and the
IDs. The general process is as follows:

1. Each __Query ID__ can "cast one vote" in total. If the Mapping
   Matrix has it voting for more than one ID (ie, it is assigned to
   two or more IDs), then the vote is evenly divided amongst the IDs.
   * For example, if probest 1234_at is assigned to four locus IDs,
     then it contributes a fractional count of 0.25 (1/4) to each of
     the four genes.
1. An ID can have at most a count of 1. That is, if mutltiple Query
   IDs are assigned to the same ID, the aggregate sum of that ID can
   not exceed 1.0
   * For example, locus XYZ3 has 3 probesets: 111_at, which
     contributes a count of 0.5 (it was also mapping to one other
     gene), 222_at provides a count of 0.25 (had mapped to three other
     genes) and 333_at providing a count of 0.33 (maps to two
     additional genes). If only 111_at and 222_at are present in a
     Query List, then XYZ3 will be tallied as having 0.75 counts
     (0.5 + 0.25). If, however, all three probesets are in a list,
     then XYZ3 is tallied as a count of 1: Although 0.25 + 0.33 + 0.75
     = 1.33, the value will be pmin()'ed to 1
1. These fractional counts are held for the entire Mapping Matrix as
   an internal matrix held in field `@mapWeights`.
    * The column sums of this matrix (`pmin()`ed to be at most 1)
      represent the maximum attainable count for any given __ID__, and
      are stored internally in `@idCount`.
1. ID counts are carried through to __Ontology Terms__, which are also
   initially expressed fractionally.
1. Finally, the fractional values are discretized before being used in
   a "normal" hypergeometric calculation. This is done with an
   internal function, `generousRound()`, which in turn relies on a
   user-configurable parameter `round`, provided as an
   integer. `round` has a default value of 2, and represents the
   denominator of the smallest fraction that will be rounded up. For
   example, a value of 4 indicates that any count `>= 1/4` should be
   rounded up, such that 17.3 will be discretized to 18, but 17.2 will
   result in 17.
    * The sum of `@idCount`, roundeded using `generousRound()`,
      provides the final discretized world size, `W`, captured as
      `@worldSize`.

For example, using HG_U133A as the Query __Namespace__ and
GeneOntology as that of the Query Terms:

* Neurofibromin 1 (NF1) has seven probesets that might be expectd to
  hybridize to its transcripts: 204323_x_at, 204325_s_at, 210631_at,
  211094_s_at, 211914_x_at, 212676_at, 212678_at
    * NF1 is only represented by those 7 probesets, and those
      probesets are each only assigned to NF1.
    * If any of those probesets are in a list, NF1 will receive a
      "full" count (1.0), as will any terms assigned to it
    * If more than one probeset is present, NF1 is still only tallied
      as 1.
* 208307_at is expected to hybridize to six distinct RBMY1 family
  members, a Y chromosome repeat family: RBMY1A1, RBMY1B, RBMY1D,
  RBMY1E, RBMY1F, RBMY1J
    * If 208307_at is present in a list, then all six loci are
      presumed "hit". This means that each ontology term assigned to
      one of these genes will "earn" 1/6th of a count from 208307_at
    * GO:0003723, "RNA binding" is assigned to all six loci. So
      208307_at will contribute a full count (6/6) toward that GO
      term.
    * GO:0007283, "spermatogenesis", is assigned to only three loci
      (B, F and J). This term will then only get a fractional count
      of 0.5 (3/6) from 208307_at
      
### Additional Features

* Multiple testing adjustment, using `p.adjust()`, defaulting to
  Benjamini & Yekutieli
* Tabular reporting using `DynamicTable`


## Scott Suggestions

* Perform deep permutation of parameter space to empirically observe
  impact on outcome.
