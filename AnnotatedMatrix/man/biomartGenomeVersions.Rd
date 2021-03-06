% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supportFunctions.R
\name{biomartGenomeVersions}
\alias{biomartGenomeVersions}
\title{Current BioMart Genome Versions}
\usage{
biomartGenomeVersions(...)
}
\arguments{
\item{\dots}{Passed to \link{biomartVersion}}
}
\value{
A vector of genome build tokens (eg "") with dataset names
    as vector names
}
\description{
Return a vector of genome versions in the current version of BioMart
}
\details{
Recovers the genomic versions used in BioMart, keyed by dataset
name associated with each.
}
\seealso{
\link{biomartVersion}
}
