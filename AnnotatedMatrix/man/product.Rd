% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodsDocumentation.R
\name{product}
\alias{product}
\title{Matrix Product}
\arguments{
\item{mat2}{Required, the second ('right') annotated matrix object}

\item{dim1}{Default \code{NULL}, the shared dimension on the left
matrix, processed by \link{sharedDimensions}.}

\item{dim2}{Default \code{NULL}. Like \code{dim1}, but for the
second (right) matrix.}

\item{ignore.case}{Default \code{TRUE}, which will intersect the
matrices without regard to case.}

\item{warn}{Default \code{TRUE}, passed to \link{sharedDimensions}}

\item{fail.na}{Default \code{TRUE}, passed to \link{sharedDimensions}}

\item{valfunc}{Default \code{NULL}, which will default to
'maxright' AND chide you to explicitly pick a value. See the
'Value Function' section for more information.}

\item{levels}{Default \code{NULL}, optional list of factor
levels. If NULL, and if valfunc is using a simple request (eg
'right max') to specify that values should come from just one
matrix, then the levels from that matrix will be
used. Otherwise you should provide a character vector of levels
if appropriate (and assure that valfunc is generating
appropriate integer output)}

\item{ignore.zero}{Default \code{TRUE}, which will report a value
of zero for any left/right pair if all values on either side
are only zero. In these situations \code{valfunc} will be
ignored.}

\item{help}{Default FALSE. If TRUE, show this help and perform no
other actions.}
}
\description{
Chain two matrices together based on a common dimension
}
\details{
\preformatted{
## Method Usage:
myObject$product( help=TRUE )

myObject$product(mat2, dim1=NULL, dim2=NULL, valfunc=NULL, levels=NULL,
                 ignore.zero=TRUE, ignore.case=TRUE)
}

Designed to chain / bridge two matrices that share a common
dimension. For example, if you have a \code{Foo vs Bar} matrix, and
a \code{Pop vs Bar} matrix, you can take their 'product' by finding
shared entities in the 'Pop' dimension, resulting in a \code{Foo vs
Pop} matrix.

Because each AnnotatedMatrix is really representing a miniature
graph database, this isn't a formally traditional matrix product,
but rather a graph traversal across the two matrices. The primary
challenge is numericly representing the final score assigned to
each Foo-vs-Pop connection, since the scores held by the matrices
are rarely going to make sense as a product. Additionally, two or
more paths will be available for any Foo-Pop combination if they
can traverse two or more intermediate Bar entries. Normalizing /
Aggregating those multiple scores into a single useful value will
not have a generic solution.
}
\section{Value Function}{


The \code{valfunc} defines a function that controls how the matrix
value will be calculated when multiple paths are available between
two output destinations. It can either be a function with two
numeric arguments reprsenting a vector of 'left' values and a
vector of 'right' values, returning a single numeric value, for
example:

\code{
   function(lv,rv) {
      max(c(lv,rv)) # Takes the maximum of either left or right value
   }
}

... or it can be a text string, which will be parsed by regular
expression:

\itemize{
  \item 'max' will return the maximum of the values
  \item 'min' returns the minimum
  \item 'mean' returns the mean
  \item 'left' or 'lft' will consider only values from the left input
  \item 'right' or 'rgt' will consider only right matrix values
  \item If neither 'right' nor 'left' (nor an alias) are detected, then both matrices values will be combined
}

So the default 'right max' will represent the value of each new
transitive pair as \code{max(rightMatrixValues)}.
}

\section{Value Function Musings}{


Consider the following two matrices; We wish to make a 'foo' ->
'bar' matrix by using the shared 'boop' dimension:

\code{
        'boop'          'boop'
         X Y             X  Y
 'foo' A 2 .     'bar' Q 5  7
       B . 3           R . 11
       C 1 4           
}

The available 'connections' are then:

\code{
    foo boop bar valRgt valLft
      A   X    Q      2      5
      B   Y    Q      3      7
      B   Y    R      3     11
      C   X    Q      1      5
      C   Y    Q      4      7
      C   Y    R      4     11
}

Partially collapsing to remove the 'discarded' boop dimension we
then would have the following distinct output dimensions, with
their [left,right] score pairs:

\code{
   A -> Q  [ [2, 5] ]
   B -> Q  [ [3, 7] ]
   B -> R  [ [3,11] ]
   C -> Q  [ [1, 5], [4, 7] ]
   C -> R  [ [4,11] ]
}

Each of these poses a challenge; We have a minimum of two numeric
values (one left, one right), and in the case of C->Q there were
two paths found, for four values total.

The matrix values could represent a wide variety of measurements;
p-values, percent similarities, rank order, factorized categories,
or simply boolean 1/0 unqualified flags. Because of this, it is
impossible to naively assure a 'relevant' mechanism to collapse /
summarize the values.

The default, 'max right', presumes that 'bigger is better' (max),
and that the second matrix is probably capturing the "final" value
we wish to use. It's a not-unreasonable approach, but given the
diversity of available approaches, will be inappropriate in many
cases.
}

\examples{

# The `machines` toy matrix associates machines with components 
mach <- AnnotatedMatrix( annotatedMatrixExampleFile("Machines.ijx") )
# The `components`  matrix associates them with parts
comp <- AnnotatedMatrix( annotatedMatrixExampleFile("Components.ijx") )

# We can use $product to associate machines with parts. These
# particular matrices have scores that are counts, so "traditional"
# matrix math is how we should set up scores:
matProd <- function(l, r) sum( l * r )

# Now we can calcualte the 'transitive' product of the two matrices:
prod <- mach$product(comp, valfunc=matProd)

# Inspect the input matrices as well as the product:
mach$matObj()
comp$matObj()
prod$matObj()

# `prod` now shows how many parts are in each machine, having
# 'joined' via the shared "Component" dimension.

}
