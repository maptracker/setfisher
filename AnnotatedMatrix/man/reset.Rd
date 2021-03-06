% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodsDocumentation.R
\name{reset}
\alias{reset}
\title{Reset Filters}
\arguments{
\item{help}{Default FALSE. If TRUE, show this help and perform no
other actions.}
}
\value{
NA, invisibly
}
\description{
AnnotatedMatrix object method to revert any applied filters
}
\details{
\preformatted{
## Method Usage:
myObject$reset( help=TRUE )

myObject$reset()
}

AnnotatedMatrix stores two matrices: \link{matrixRaw}, which
remains an unaltered form of the matrix as loaded from the source
file, and \link{matrixUse}, the version resulting from application
of filters. The reset() method will set matrixUse to NULL, and
clear the $setFilters and $filterLog structures.
}
\examples{

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )

## The example matrix has a default filter that suppresses a few symbols:
sketchySymbol <- "B(p51A)"
# No output, though Score=0 (not NA) shows symbol is recognized
s2e$map( sketchySymbol )  

s2e$reset()
s2e$map( sketchySymbol )  # Connection reestablished to LOC8626

}
\seealso{
\link{matObj}, \link{matrixRaw}, \link{matrixUse},
    \link{filterLog}, \link{filterSummary}, \link{appliedFilters},
    \link{autoFilter}
}
