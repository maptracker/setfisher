## RefClass objects have very minimal support for documentation in
## ROxygen (possibly little support at all?). This document is an
## attempt to centralize method documentation. Preferably these blocks
## would be co-linear with the actual function, but RefClass
## $methods() specification makes that difficult / impossible for all
## but the first method.


#' Raw Matrix
#'
#' Internal AnnotatedMatrix field holding unaltered (as read from
#' file) Matrix object
#'
#' @name matrixRaw
#'
#' @details
#'
#' This AnnotatedMatrix object field is populated directly with matrix
#' data from the source file, and is unaltered. It allows filters to
#' be cleared using \link{reset}.
#'
#' @return a Matrix:: 'dgTMatrix' sparse matrix
#'
#' @seealso \link{matObj}, \link{matrixUse}, \link{reset}
#'
#' 
NULL

#' Used (Filtered) Matrix
#'
#' Internal AnnotatedMatrix field holding the filtered Matrix object
#'
#' @name matrixUse
#'
#' @details
#'
#' This AnnotatedMatrix object field is \code{NULL} if no filters have
#' been applied. Otherwise it is a filtered form of
#' \link{matrixRaw}. You will generally not want to access this field
#' directly; Instead, use \link{matObj}, as it will by default provide
#' matrixUse if filters are applied, or will fallback to
#' matrixRaw. Using \link{reset} will reset all filters by setting
#' matrixUse to \code{NULL}.
#'
#' @return a Matrix:: 'dgTMatrix' sparse matrix, or NULL if no filters
#'     have been applied
#'
#' @seealso \link{matObj}, \link{matrixRaw}, \link{reset}
#'
#' 
NULL

#' Filter Log
#'
#' Internal AnnotatedMatrix field holding a data.table of rows or
#' columns removed by filters
#'
#' @name filterLog
#'
#' @details
#'
#' This AnnotatedMatrix object field is \code{NULL} if no filters have
#' been applied. Otherwise it is a filtered form of
#' \link{matrixRaw}. You will generally not want to access this field
#' directly; Instead, use \link{matObj}, as it will by default provide
#' matrixUse if filters are applied, or will fallback to
#' matrixRaw. Using \link{reset} will reset all filters by setting
#' matrixUse to \code{NULL}.
#'
#' @return A data.table
#'
#' @seealso \link{reset}, \link{filterSummary}
#'
#' @examples
#'
#' # The example matrix includes some automatic filters:
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$filterSummary()  # Overview of filter effects
#' # What columns have been eliminated?
#' s2e$filterLog[ s2e$filterLog$type == 'Col' , c("id", "metric") ]
#'
#' # Frank says he can't find data for p51 anymore. Was it filtered?
#' s2e$filterLog[ s2e$filterLog$id == "p51", ]
#' # Yes; Provide Frank with guidance on avoiding sketchy gene symbols
#' 
NULL

#' Get Matrix Object
#'
#' AnnotatedMatrix object method to recover the 'current' matrix data
#'
#' @name matObj
#' @method matObj AnnotatedMatrix
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myAnnMat$matObj( help=TRUE )
#' 
#' obj <- myAnnMat$matObj( raw=FALSE, transpose=FALSE )
#' }
#'
#' @param raw Default \code{FALSE}, which will return \link{matrixUse}
#'     if filters have been applied. If no filters are applied, or
#'     raw=TRUE, then \link{matrixRaw} will be returned
#' @param transpose Default \code{FALSE}. if TRUE, then transpose the
#'     matrix before returning it.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return a Matrix:: 'dgTMatrix' sparse matrix
#'
#' @seealso \link{reset}, \link{matrixRaw}, \link{matrixUse}
#'
#' @examples
#' 
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' ## Automatically applied filters result in fewer non-zero elements
#' ## of the filtered matrix than the raw one:
#' Matrix::nnzero(s2e$matObj())         # 180
#' Matrix::nnzero(s2e$matObj(raw=TRUE)) # 208
#'
#' # transpose is included to help you avoid using base::t(), which
#' # will de-sparsify the matrix and often result in RAM exhaustion
#' head( rownames( s2e$matObj() ) )
#' head( rownames( s2e$matObj( transpose=TRUE ) ) )
#' 
NULL

#' Row Names
#'
#' AnnotatedMatrix object method to get/set/reorder the rownames of
#' the filtered matrix
#'
#' @aliases rNames cNames
#' @method rNames AnnotatedMatrix
#' @method cNames AnnotatedMatrix
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myAnnMat$rNames( help=TRUE )
#' 
#' rn <- myAnnMat$rNames( raw=FALSE )
#' cn <- myAnnMat$cNames( raw=FALSE )
#'
#' myAnnMat$rNames( newRowNames, reason=NA )
#' myAnnMat$cNames( newColNames, reason=NA )
#' }
#'
#' Normally would be used to simply recover the rownames for the
#' active matrix. Can also be used to reorder, or to reduce the
#' accessible names. Unrecognized names will be added to the matrix as
#' a row or column of all zeros; In most cases this will not be useful
#' (no connection to other dimension) but might be helpful in some
#' circumstances.
#'
#' @param new Default NULL, which will simply return the current row
#'     or column names. Alternatively can provide a character vector
#'     of new values. This can be used to reorder, remove or add
#'     names. Added names will generate a row or column of all zeros
#'     (no connections). Setting this value is treated as a filter
#'     action, so \link{matrixUse} will be altered and actions noted
#'     in \link{filterLog}.
#' @param raw Default FALSE, which will operate on the filtered
#'     matrix. See \link{matObj} for more information. A fatal error
#'     will be thrown if raw is TRUE and new is not NULL (can not
#'     alter raw matrix)
#' @param reason Default NA. Only relevant when row is not
#'     NULL. Optional human-readable reason for the alteration was
#'     made, will be recorded in \link{filterLog} and used to
#'     structure \link{filterSummary}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character vector of row names
#'
#' @seealso \link{reset}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' length( s2e$rNames() ) # 167 names known
#' s2e$map('LOC326')      # Maps to six symbols
#' s2e$rNames(c('AIRE1', 'PGA1', 'AngryPirate', 'AIRE'),
#'            reason="Nancy asked me to limit analysis to these four genes")
#' s2e$filterSummary()    # Alteration is recorded in filterLog
#' length( s2e$rNames() ) # We now have the above four
#' s2e$map('LOC326')      # We only kept three of the original matches
#' s2e$removeEmpty()      # Remove rows (and columns) that are only zeros
#' s2e$matObj()           # A much smaller matrix!
#'
NULL

#' Reset Filters
#'
#' AnnotatedMatrix object method to revert any applied filters
#' 
#' @name reset
#' @method reset AnnotatedMatrix
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myAnnMat$reset( help=TRUE )
#' 
#' myAnnMat$reset()
#' }
#'
#' AnnotatedMatrix stores two matrices: \link{matrixRaw}, which
#' remains an unaltered form of the matrix as loaded from the source
#' file, and \link{matrixUse}, the version resulting from application
#' of filters. The reset() method will set matrixUse to NULL, and
#' clear the $setFilters and $filterLog structures.
#'
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return NA, invisibly
#'
#' @seealso \link{matObj}, \link{matrixRaw}, \link{matrixUse},
#'     \link{filterLog}, \link{filterSummary}, \link{appliedFilters}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' 
#' ## The example matrix has a default filter that suppresses a few symbols:
#' sketchySymbol <- "B(p51A)"
#' # No output, though Score=0 (not NA) shows symbol is recognized
#' s2e$map( sketchySymbol )  
#'
#' s2e$reset()
#' s2e$map( sketchySymbol )  # Connection reestablished to LOC8626
#'
NULL

#' Automatic Filters
#'
#' AnnotatedMatrix object method to apply filters described in parameters
#'
#' @name autoFilter
#' @method autoFilter AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myAnnMat$autoFilter( help=TRUE )
#'
#' myAnnMat$autoFilter( recursive=TRUE, verbose=TRUE )
#' }
#' 
#' The matrix will recognize several parameters as defining filter
#' operations. While these parameters could be set 'manually' in code,
#' they are of primary utility in allowing annotated matrix files to
#' "self describe" default reasonable starting filters for the matrix.
#'
#' The following parameters are identified as specifying a filter
#' request:
#' 
#'   MinScore   MaxScore
#'   KeepLevel  TossLevel
#'   TossMeta
#'
#' Calling \code{$showParameters( na.rm=FALSE)} will list all standard
#' parameters and their values; Not all of these will be relevant to
#' autoFilter, but those that are will indicate as much in the
#' description.
#'
#' @param recursive Default TRUE. Because some filters check the
#'     number of populated rows or columns, it is possible that the
#'     application of a latter filter will result in rows or columns
#'     that passed a prior one to now be considered failing. Recursive
#'     will keep reapplying filters into no further changes (based on
#'     zeroed cell count) are observed. If you desire the final matrix
#'     to pass all filters, leave recursive as TRUE.
#' @param verbose Default TRUE, which displays a message summarizing
#'     the count of zeroed elements (cells, rows and columns).
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return Invisibly, an integer vector of newly zeroed counts for
#'     (cells, rows, cols).
#'
#' @seealso \link{filterByScore}, \link{filterByFactorLevel},
#'     \link{filterByMetadata}, \link{reset}, \link{filterSummary},
#'     \link{appliedFilters}, \link{showParameters}
#'
#' @examples
#'
#' ## The example matrix has two default filters defined:
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e # Summary of the matrix
#' 
#' ## Only one removed any entries, though:
#' s2e$filterSummary()
#' ## This is because the cells impacted by the second filter (on
#' ## "Description" metadata) are fully contained by those already
#' ## removed by the first (a factor filter removing level "Unknown")
#'
#' ## We can remove all filters:
#' s2e$reset()
#' s2e
#'
#' ## ... and then reapply the automatic filters if we wish:
#' s2e$autoFilter()
#'
#' ## Show the currently set parameters:
#' s2e$showParameters()
#'
#' ## Show all recognized parameters:
#' s2e$showParameters( na.rm=FALSE )
#'
NULL
