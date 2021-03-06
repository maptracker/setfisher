% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodsDocumentation.R
\name{DOTbuildMatrix}
\alias{DOTbuildMatrix}
\alias{.buildMatrix}
\title{Build Matrix}
\arguments{
\item{metadata}{Default \code{NULL}. Optional metadata
table. Should be a \code{data.table} keyed to column 'id'.}

\item{levels}{Default \code{NULL}. Optional character vector of
levels. Will cause the matrix to be treated as a pseudo-factor.}

\item{cols}{Default \code{NULL}. Optional metadata column
definitions. Should be a list with column names as keys,
human-readable descriptions as values.}

\item{help}{Default FALSE. If TRUE, show this help and perform no
other actions.}
}
\value{
At the moment, just \code{TRUE}.
}
\description{
Internal AnnotatedMatrix object method to generate matrix 'by hand'
}
\details{
Internal method, should not be called directly. It is called
automatically when a new object is created with
\code{AnnotatedMatrix(someDataStructure)}.
}
\seealso{
\link{.readMatrix}, \link{.readObject}
}
