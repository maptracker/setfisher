## RefClass objects have very minimal support for documentation in
## ROxygen (possibly little support at all?). This document is an
## attempt to centralize method documentation. Preferably these blocks
## would be co-linear with the actual function, but RefClass
## $methods() specification makes that difficult / impossible for all
## but the first method.


#' Raw Matrix
#'
#' Internal field holding unaltered (as read from file) Matrix object
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
#' Internal field holding the filtered Matrix object
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

#' Get Matrix Object
#'
#' Recover the 'current' matrix data
#'
#' @name matObj
#'
#' @param raw Default \code{FALSE}, which will return \link{matrixUse}
#'     if filters have been applied. If no filters are applied, or
#'     raw=TRUE, then \link{matrixRaw} will be returned
#' @param transpose Default \code{FALSE}. if TRUE, then transpose the
#'     matrix before returning it.
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
NULL

#' Reset Filters
#'
#' Revert any applied filters such that the 'used' matrix is same as 'raw'
#' 
#' @name reset
#'
#' @details
#'
#' AnnotatedMatrix stores two matrices: \link{matrixRaw}
#'
#' @return NA, invisibly
#'
#' @seealso \link{matObj}, \link{matrixRaw}, \link{matrixUse}
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
