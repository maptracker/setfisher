## ROxygen provides only rudimentary support for documenting Reference
## Class object methods (as of 2017). It will recognize an unassigned
## character string at the beginning of the method function as a
## description, and will parse the parameters. This class is complex
## enough that I wanted more formalized documentation. The blocks
## below define help topics for both RefClass fields and methods:

## * The method name is being put in @name and/or @aliases. This will
##   likely be interpreted by the help system as a simple function in the
##   namespace, but I don't think that will cause problems. It allows
##   \link{}s to be built between topics.
## * I am also specifying the methods under the @method section. This
##   will likely be perceived as an S3/S4 function in the class. I don't
##   _think_ that's an issue, but it might be.
## * It's effectively impossible to use the @usage section, since R_CMD
##   check becomes very unhappy with attempts to formalize usage in
##   'true' object fasion. Instead, usage is manually smuggled into the
##   @details section in a \preformatted{} block.

## ======================= FIELDS =======================

#' Default Query Object
#'
#' Internal SetFisher field holding the Default Query Matrix
#'
#' @name defQryObj
#'
#' @details
#'
#' Either NULL (technically an empty field object) if no default has
#' been defined, or an AnnotatedMatrix object. Set using
#' \link{defaultQuery}.
#'
#' @seealso \link{defaultQuery}
#'
#' @examples
#'
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultQuery( s2e )
#' str( sf$defQryObj )
NULL

#' Default Map Object
#'
#' Internal SetFisher field holding the Default Map Matrix
#'
#' @name defMapObj
#'
#' @details
#'
#' Either NULL (technically an empty field object) if no default has
#' been defined, or an AnnotatedMatrix object. Set using
#' \link{defaultMap}.
#'
#' @seealso \link{defaultMap}
#'
#' @examples
#'
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultMap( s2e )
#' str( sf$defMapObj )
NULL

#' Default Query World
#'
#' Internal SetFisher field holding the default set of query IDs
#'
#' @name defQueryWorld
#'
#' @details
#'
#' Either NULL (technically an empty field object) if no default has
#' been defined, or a character vector of Query IDs. Set with
#' \link{defaultQueryWorld}.
#'
#' @seealso \link{defaultQueryWorld}, \link{defaultQuery}
#'
#' @examples
#'
#' sf  <- SetFisher()
#' world <- c("ABC", "XYZ", "ZF12", "HOX4", "OCT4")
#' sf$defaultQueryWorld( world )
#' str( sf$defQueryWorld )
NULL



## ======================= METHODS ======================

#' Default Query Matrix
#'
#' Set/Get an optional default Query Matrix
#'
#' @name defaultQuery
#' @method defaultQuery SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$defaultQuery( help=TRUE )
#' 
#' queryMat <- myObject$defaultQuery( x=NULL, ... )
#' }
#'
#' If a default query matrix is set, it will be automatically used in
#' cases where one is needed but none is supplied (generally when
#' creating a \link{SetFisherAnalysis})
#'
#' The query is a set of one or more sets of IDs that will be tested
#' for enrichment in one or more ontologies.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     query matrix, or NULL if none is set. Otherwise, will attempt
#'     to recover the AnnotatedMatrix corresponding to \code{x} and
#'     set it as the default.
#' @param \dots Passed to \link{annotatedMatrix}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NULL} if none set, or an \link{AnnotatedMatrix}
#'
#' @seealso \link{defQryObj}, \link{defaultQueryWorld}, \link{defaultMap},
#'     \link{annotatedMatrix}, \link{AnnotatedMatrix}
#'
#' @importFrom CatMisc is.empty.field
#'
#' @examples
#' 
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultQuery( s2e )
#' 
NULL

#' Default Map Matrix
#'
#' Set/Get an optional default Map Matrix
#'
#' @name defaultMap
#' @method defaultMap SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$defaultMap( help=TRUE )
#' 
#' mapMat <- myObject$defaultMap( x=NULL, ... )
#' }
#' 
#' If a default mapping matrix is set, it will be automatically used
#' in cases where one is needed but none is supplied (generally when
#' creating a \link{SetFisherAnalysis})
#' 
#' The map is an optional matrix that converts IDs from one
#' "namespace" to another. It is used to connect queries to ontologies
#' when they do not share a namespace. For example, if your queries
#' are sets of Affymetrix probe sets, but your ontology is defined
#' using Ensembl gene accessions, you will need a map that associates
#' Affy probe set IDs with the appropriate Ensembl Gene IDs.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     map matrix, or NULL if none is set. Otherwise, will attempt
#'     to recover the AnnotatedMatrix corresponding to \code{x} and
#'     set it as the default.
#' @param \dots Passed to \link{annotatedMatrix}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NULL} if none set, or an \link{AnnotatedMatrix}
#'
#' @seealso \link{defMapObj}, \link{defaultQuery},
#'     \link{annotatedMatrix}, \link{AnnotatedMatrix}
#'
#' @importFrom CatMisc is.empty.field
#'
#' @examples
#' 
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultMap( s2e )
#' 
NULL

#' Default Query World
#'
#' Set/Get an optional default Query World
#'
#' @name defaultQueryWorld
#' @method defaultQueryWorld SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$defaultQueryWorld( help=TRUE )
#' 
#' myWorld <- myObject$defaultQueryWorld( x=NULL )
#' }
#' 
#' If a default query world is set, it will be automatically used
#' in cases where one is needed but none is supplied (generally when
#' creating a \link{SetFisherAnalysis})
#'
#' The query world is used to constrict or expand the collection of
#' \strong{all possible} query IDs. For example, you may be running an
#' analysis on mouse genes, but your experimental proceedure
#' explicitly sampled a smaller subset of specific genes. In this
#' case, it is important to constrain the world to that set of genes
#' to properly calculate enrichment statistics.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     query world, or NULL if none is set. Otherwise, will expect a
#'     character vector of IDs that defines the world of query IDs you
#'     wish to use.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NA} if none set, or a character vector of IDs that
#'     defines the world being analyzed.
#'
#' @seealso \link{defQryWorld}, \link{defaultQuery}
#'
#' @importFrom CatMisc is.empty.field
#'
#' @examples
#' 
#' sf  <- SetFisher()
#' world <- c("ABC", "XYZ", "ZF12", "HOX4", "OCT4")
#' sf$defaultQueryWorld( world )
NULL

#' Get AnnotatedMatrix
#'
#' Get an AnnotatedMatrix, either from a cache or de novo
#'
#' @name annotatedMatrix
#' @method annotatedMatrix SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$annotatedMatrix( help=TRUE )
#' 
#' myMat <- myObject$annotatedMatrix( x=NA, ... )
#' }
#' 
#' This is a light-weight wrapper around an \code{AnnotatedMatrix()}
#' call, with some additional functionality that an already-loaded
#' matrix will be recovered from an internal cache rather than
#' reloaded.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     query world, or NULL if none is set. Otherwise, will expect a
#'     character vector of IDs that defines the world of query IDs you
#'     wish to use.
#' @param \dots Passed to \link{AnnotatedMatrix}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NULL} if none set, or a character vector of IDs that
#'     defines the world being analyzed.
#'
#' @seealso \link{AnnotatedMatrix}, \link{matrixCache}
#'
#' @import AnnotatedMatrix
#'
#' @examples
#'
#' \dontrun{
#' sf  <- SetFisher()
#' myMat <- sf$annotatedMatix("/data/matrices/Chickens_vs_Cows.mtx")
#' }
#' 
NULL

