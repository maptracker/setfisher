% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supportFunctions.R
\name{DOTmakeTempFile}
\alias{DOTmakeTempFile}
\alias{.makeTempFile}
\title{Make Temp File}
\usage{
.makeTempFile(name, pkg = "AnnotatedMatrix")
}
\arguments{
\item{name}{Required, the basename of the file}

\item{pkg}{Default "AnnotatedMatrix", the package name}
}
\description{
Internal method to help manage example files in extdata/
}
\details{
When matrix files are loaded, a .rds version is generated to aid in
rapid loading in the future. These files are made alongside the
original flat file. This could cause problems if the file is from
inst/extdata, so this method takes a requested file and copies it
to a temporary location before loading. It will also copy sidecar
files used by the 'primary' file.
}
\examples{

tmpFile <- AnnotatedMatrix:::.makeTempFile("Symbol-To-Gene.mtx")

}
\seealso{
\link{sidecarFiles}
}
