SetFisher is an R package that manages the application of `phyper()`
for [hypergeometric distribution][HD] ([Fisher's exact test][FET])
analyses supporting [gene set enrichment analysis][GSEA]. Novel
components of the package are:

* __Namespace mapping__: An optional "Mapping Matrix" can be provided,
  which maps input IDs to different identifiers used in the ontologies
  (eg from Affymetrix probesets to Entrez GeneIDs)
* __Multiple voting accomodation__: Fisher's Exact Test is very
  sensitive to dependencies amongst the input (it's why it is used,
  after all). If there are 1:many or many:many relationships between
  the experimental assays and the target genes, this technical
  dependency can obscure biological relationships. SetFisher attempts
  to manage these issues by tracking "fractional counts" through the
  mapping matrix.
* __Automated filtering__: Ontologies can be optionally "trimmed" to
  remove terms with few assigned genes, or those with a large
  number. Similarly, genes can be optionally filtered to remove those
  with insufficient ontological support. This process helps remove
  "speculative" genes that are either fictional (do not actually
  exist) or extremely esoteric; Such genes tend to inflate
  significance scores.

More details can be found in [How it Works][HiW], which lays out the
philosophy and basic operation of the package.

Sample analyses and ontology matrices will be published at a latter date.


[HD]: https://en.wikipedia.org/wiki/Hypergeometric_distribution
[FET]: https://en.wikipedia.org/wiki/Fisher%27s_exact_test
[GSEA]: https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis
[HiW]: HowItWorks.md
