% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodsDocumentation.R
\name{modState}
\alias{modState}
\title{Modification State}
\value{
A numeric value that is increased with each modification
}
\description{
Internal AnnotatedMatrix field tracking modification events
}
\details{
The "raw" matrix (\link{matrixRaw}) will always stay
unaltered. However some operations, primarily filters, will change
the "used" matrix (\link{matrixUse}). It may be of value for
external code to recognize when an AnnotatedMatrix has undergone
such an alteration. The \code{$modState} field is designed to allow
that; If the value held by external code differs from what the
matrix reports, then the matrix has since undergone an operation
that will (almost certainly) have changed it (in some unspecified
way).

\preformatted{
## THIS FIELD SHOULD NOT BE MODIFIED BY THE USER OR OTHER CODE
## DOING SO COULD LEAD TO CONFUSION / CRASHES
}
}
\examples{

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
s2e$modState # The initial state. Autofilters may have incremented this.

s2e$cNames()
s2e$modState # no change, we just looked at some info

s2e$filterByFactorLevel("Official")
s2e$modState # We altered the matrix

s2e$reset()
s2e$modState # Again, altered the matrix
}
\seealso{
\link{matrixUse}
}
