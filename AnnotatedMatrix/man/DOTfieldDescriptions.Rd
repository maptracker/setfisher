% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodsDocumentation.R
\name{DOTfieldDescriptions}
\alias{DOTfieldDescriptions}
\alias{.fieldDescriptions}
\title{Field Descriptions}
\arguments{
\item{update}{Default TRUE, which will modify the internal object
field to have the two attached objects.}

\item{help}{Default FALSE. If TRUE, show this help and perform no
other actions.}
}
\value{
A list, with field variable names as names, and descriptive
    text as values.
}
\description{
Internal AnnotatedMatrix object method providing descriptive text
to object fields
}
\details{
Internal method, should not be called directly

Helper method that adds descriptive attributes to the fields held
by the object. This is an attempt to make the object
self-documenting. Two fields are added: "Description", which
provides a one-sentance description of the field, and "Help", which
is a copy-and-pastable \link{help} command to recover the man page
for the field
}
\seealso{
\link{help}
}
