% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodsDocumentation.R
\name{filterByMetadata}
\alias{filterByMetadata}
\title{Filter by Metadata}
\arguments{
\item{key}{Required, the name of the key (tag / column) to match
against. Will throw an error if not a single value}

\item{val}{Required, the value to match against. Can be multiple
values.}

\item{MARGIN}{Default NULL, which will look for matches in both
rows in columns. Alternatively can be 1, 2, "row" or "col".}

\item{keep}{Default TRUE, which will keep any matching entry,
discarding (zeroing out) the others. If FALSE, then rows/cols
with metadata which does not match will be removed.}

\item{type}{Default 'like', controls the kind of matching. 'like'
will allow a fragmentary match of your value to metadata (both
"left" and "right" wildcards), 'regexp' will treat your input
as a regular expression, and "equal" will require an exact
match. Additionaly, the text "case" can be added to any of the
types to force a case-sensitive match; For example,
"equal case" will require an exact, case-sensitive match.}

\item{op}{Default 'or', controls how filters resolve two or more
\code{val} entries. 'or' will count a metadata entry as
matching if any of your \code{val} entries match, while 'and'
requires all to match. Additionally, you can use aliases 'any'
(='or') or 'all (='and').}

\item{filterEmpty}{Default FALSE. If true, then rows or columns
that are empty after filtering will be removed from the matrix
using \link{removeEmpty}. If no cells were affected by the
\code{filterByMetadata} call, then \code{removeEmpty} will not
be called. If you wish to assure that empty rows and columns
are removed after this filter, you should call
\code{removeEmpty} explicitly.}

\item{reason}{Default NA. Optional human-readable reason for why
the alteration was made, will be recorded in \link{filterLog}
and used to structure \link{filterSummary}}

\item{help}{Default FALSE. If TRUE, show this help and perform no
other actions.}
}
\value{
Invisibly, an integer vector of newly zeroed counts for
    (cells, rows, cols).
}
\description{
AnnotatedMatrix object method to filter rows or cells by metadata
}
\details{
\preformatted{
## Method Usage:
myObject$filterByMetadata( help=TRUE )

myObject$filterByMetadata(key, val, MARGIN=NULL, type="like", op='or',
                          reason=NA)
}

Updates the \link[=matrixUse]{used matrix} to remove (zero-out) all
rows or columns that have metadata failing to match your
request. Like all filters, will not alter \link{matrixRaw}.

Note that this is not a cell-level filter, since metadata is
associated with the dimensions of the matrix, not the "connections"
between dimensions. It will thus remove entire rows or
columns. However, it will still tally the number of newly-zeroed
cells generated by the filter operation.
}
\examples{

# The example matrix has some deprecated genes. These are already
# removed by the default auto filters, so we will reset the matrix
# before looking for them
s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
s2e$reset()

# This tag is used in the Description field to highlight deprecated loci:
depTag <- "{Deprecated}"
s2e$filterByMetadata("description", depTag, keep=TRUE, filterEmpty=TRUE,
                        reason="Inspecting deprecated loci")
s2e$filterSummary()

# We set filterEmpty to TRUE in order to more easily inspect the
# matrix iteslf:
s2e$matrixUse

# Verify that we found what we were looking for:
s2e$metadata( colnames(s2e$matrixUse), "Description" )

}
\seealso{
\link{filterLog}, \link{filterSummary},
    \link{appliedFilters}, \link{metadata}, \link{autoFilter}
}
