### Template file

This is a template file that is copied into the `byNamespace`
directories. This text (above the horizontal rule) will not be copied.

----
# AnnotatedMatrix Files - Namespace Hierarchy

### CAUTION CAUTION CAUTION

Files may be loaded directly from this folder tree, but their location
is not persistent! The tree's primary use is to aid in exploration of
available files. If you are encoding a file in a production workflow,
please reference instead the [../byAuthority][BA] file paths. Those
locations are static and will persist as new versions are loaded.

### Structure

This folder contains a file structure that organizes
[MatrixMarket][MTX] files according to their row and column
namespaces. At the level of this file are folders for namespaces on
_either_ the rows or columns. The second level will be the namespaces
of the other dimension of the matrix. Within that second level will be
symlinks to the files themselves. The symlinks will point back up to
[../byAuthority][BA], which is the primary storage location for the
files. For example:

```
byNamespace/
  EntrezGene/
    EnsemblGene/
      Map@EntrezGene_to_EnsemblGene@HGNC@2017-07-26.mtx ->
         ../../../byAuthority/HGNC/2017-07-26/Map@EntrezGene_to_EnsemblGene@HGNC@2017-07-26.mtx
    MGD/
      Map@MGD_to_EntrezGene@HGNC@2017-07-26.mtx ->
         ../../../byAuthority/HGNC/2017-07-26/Map@MGD_to_EntrezGene@HGNC@2017-07-26.mtx
```

The files can be read and manipulated by the [AnnotatedMatrix][AM] R
package. The structure is designed to aid in human discovery of useful
matrices while browsing the file system. Because _both_ row and column
namespaces will be present at _each_ of the two levels, each file will
be represented twice within this structure.

### What's a Namespace?!?

A 'namespace' is a collection of identifiers; Examples include things
like zip codes, social security numbers, stock exchange ticker
symbols, Wikipedia topic names. Biological examples are Ensembl
Protein IDs, Entrez Gene IDs, gene symbols, GeneOntology terms, PubMed
IDs. Identifiers are often segregated by namespaces to prevent
'collisions' between IDs that are the same but refer to different
entities. For example, the BrainArray probesets based on Entrez Gene
IDs have many collisions with Affymetrix probesets of unrelated genes:

```
Affymetrix : 859_at = cytochrome P450 family 1 subfamily B member 1
BrainArray : 859_at = Caveolin 3
```

In addition to the `.mtx` flat files, the directory may also contain
`.mtx.rds` files. These are serialized R data structures created after
the `.mtx` file is first parsed. In general, you should continue to
work directly with the `.mtx` files - the code can detect updates
(based on date stamps) and refresh the RDS files as needed.


[MTX]: http://math.nist.gov/MatrixMarket/formats.html
[AM]: https://github.com/maptracker/setfisher/tree/master/AnnotatedMatrix
[BA]: ../byAuthority/
