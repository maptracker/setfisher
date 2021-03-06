% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supportFunctions.R
\name{inst.path.package}
\alias{inst.path.package}
\title{Development Package Path}
\usage{
inst.path.package(pkg = "AnnotatedMatrix")
}
\arguments{
\item{pkg}{Default "AnnotatedMatrix", the package name}
}
\value{
The path to the package folder, or the inst/ subfolder if
    it exists
}
\description{
Get the path to package file, including when the package is locally
installed from a local installation
}
\details{
When developing with the uncompressed R package, some files are in
different locations. In particular, exdata/ remains inside the
inst/ folder. This method detects the presence of "inst" and
includes it in the return path
}
\seealso{
\link[base]{path.package}
}
