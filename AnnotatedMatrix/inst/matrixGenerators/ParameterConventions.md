# Parameter Conventions

Parameters are an important part of the "annotation" in
AnnotatedMatrix. They provide information regarding the matrix as a
whole, and are accessed via the `$param()` object method (row and
column annotation are accessed via the `$metadata()` object method).

There are a handful of parameters that are "reserved" and are
explicitly referenced by internal code. While these are certainly
important, other parameters are also valuable in helping to describe a
matrix to a human user. Keeping some level of consistency for these
additional parameters will make it _much_ easier for software to
consistently consume and interpret the parameters.

Notes: Parameter names are case insensitive, but the system will
'remember' the case that's first used for pretty-printing. Each
parameter can also be associated with a comment, which will be shown
by the default `$show()` method to help inform the user of their
meaning.

## Reserved Parameters

These parameters are hardcoded within AnnotatedMatrix. While it is
unlikely that the system will "break" if different names are used,
doing so would prevent "expected" information from being properly
presented.

* `Name` - A 'nice' short name to label the entire matrix
* `Description` - A longer description of what the matrix is
* `ScoreDesc` - A description of what the numeric score values
  represent.
* `Source` - Intended to be a URL or other identifier to help the user
  find the primary source of the data in the matrix.
* `Authority` - The entity (person, organization, software) that
  generated the primary data held by the matrix (not the code that
  formatted the matrix)
* `Version` - A version ID for the primary data in the matrix
* `RowDim` - The name for the dimension represented by row identifiers
* `ColDim` - The name for the dimension represented by column
  identifiers
* `RowUrl` - A sprintf format used to build URLs for specific Row IDs
* `ColUrl` - A sprintf format used to build URLs for specific Col IDs
* `CellMetadata` - If present, specifies that the cell values are
  object IDs that have additional metadata associated with them.
* `CellDim` - The name for the dimension represented by cell
  identifiers
* `LoadComment` - A message to show when the matrix is loaded.

### Reserved AutoFilter Parameters

These parameters will be interpreted as default filter settings 

* `MinScore`, `MaxScore` - Filter cells based on their numeric score
* `KeepLevel`, `TossLevel` - Filter cells based on their factor level
  (non-sensical if the matrix is not factorized)
* `KeepRow`, `TossRow`, `KeepCol`, `TossCol` - Remove (or keep only)
  rows or cols with specific IDs
* `MinRowCount`, `MaxRowCount`, `MinColCount`, `MaxColCount` - Remove
  (or keep only) rows or cols with more or less than a threshold
  number of assigned cells.
* `AutoFilterComment` - Message to show when `$autoFilter()` runs (eg
  to explain logic used, or just remind user that some data has been
  automatically suppressed)

## Other Parameters

AnnotatedMatrix will effectively ignore all other parameters. However,
it is strongly encouraged that you adopt some conventions for the
names you use and stick with them, to aid in interoperability between
matrices. Here are some of the other parameters we are using in our work:

* `Citation` - A citation for the primary source, usually specified by
  the authority that provides the data.
* `Revision` - A version number for matrix generation. That is,
  `Version` is used to reflect the state of the primary data, while
  `Revision` is used to reflect changes in the code that process the
  primary data into an AnnotatedMatrix.
* `Modifier` - Kind of like a dimension for the whole matrix. For
  example, a matrix that converts human Entrez Gene to human Ensembl
  Gene would have a RowDim of "Entrez Gene", a ColDim of "Ensembl
  Gene", and a Modifier of "Human". Alternatively, there could be a
  human Entrez Gene to mouse Entrez Gene matrix. Then RowDim="Human",
  ColDim="House mouse", and Modifier="Entrez Gene"
* `Species` - The scientific name for the species associated with the
  matrix, eg "Homo sapiens" or "Mus musculus".
