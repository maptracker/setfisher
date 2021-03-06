% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supportFunctions.R
\name{biomartCommonNames}
\alias{biomartCommonNames}
\title{BioMart Common Names}
\usage{
biomartCommonNames(...)
}
\arguments{
\item{\dots}{Passed to \link{biomartVersion}}
}
\value{
A vector of species common names (eg "Crab-eating macaque")
    with dataset names as vector names
}
\description{
Dataset-to-Common-Name lookup, providing "human friendly" species names
}
\details{
Recovers the species common names used in BioMart, keyed by dataset
name associated with each.
}
