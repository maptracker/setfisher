## ROxygen provides only rudimentary support for documenting Reference
## Class object methods (as of 2017). It will recognized an unassigned
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

#' Matrix File Path
#'
#' Internal AnnotatedMatrix field holding the path to the matrix file
#'
#' @name file
#'
#' @details
#'
#' This is simply the path to the raw file (or directory of files)
#' that define the matrix.
#'
#' @return A character vector with a single value
#'
#' @seealso \link{fromRDS}, \link{.readMatrix}
#'
#' @examples
#' 
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$file
#' 
NULL

#' Matrix Source RDS Flag
#'
#' Internal AnnotatedMatrix field defining if loading was from a cached RDS
#'
#' @name fromRDS
#'
#' @details
#'
#' A logical flag, if TRUE it indicates that the matrix was loaded
#' from a cached RDS object. The RDS file will always be the same as
#' \link{file}, but with an additional '.rds' suffix appended.
#'
#' If the RDS file does not already exist, \link{.readMatrix} will
#' parse the 'original' files and automatically create a new RDS
#' file. Similarly, if the RDS file is older than the original
#' sources, it will be regenerated from source.
#'
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' }
#'
#' @return A character vector with a single value
#'
#' @seealso \link{file}, \link{.readMatrix}
#'
#' @examples
#' 
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$fromRDS
#' 
NULL

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
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' ## Use methods instead, in this case $matObj()
#' ## ALTERING THIS FIELD WILL LEAD TO CHAOS. PLEASE DON'T.
#' }
#'
#' @return a Matrix:: 'dgTMatrix' sparse matrix
#'
#' @seealso \link{matObj}, \link{matrixUse}, \link{reset}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' str(s2e$matrixRaw)
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
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' ## Use methods instead, in this case $matObj()
#' ## ALTERING THIS FIELD WILL LEAD TO CHAOS. PLEASE DON'T.
#' }
#'
#' @return a Matrix:: 'dgTMatrix' sparse matrix, or NULL if no filters
#'     have been applied
#'
#' @seealso \link{matObj}, \link{matrixRaw}, \link{reset}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' # Example matrix has an autofilter, so matrixUse is populated:
#' str(s2e$matrixUse)
#' # Resetting filters will also clear $matrixUse, until new filters
#' # are applied.
#' s2e$reset()
#' str(s2e$matrixUse)
#' 
NULL

#' Matrix Metadata Object
#'
#' Internal AnnotatedMatrix field holding metadata about row and column entries
#'
#' @name matrixMD
#'
#' @details
#'
#' This field stores metadata information for the AnnotatedMatrix
#' object. Metadata for rows and columns are both stored in a single
#' \code{data.table} (this may be a bad idea - I may split the field
#' later into matrixMdRow and matrixMdCol)
#'
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' ## Use methods instead, in this case $metadata()
#' }
#' @return a data.table object
#'
#' @seealso \link{metadata}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' str(s2e$matrixMD)
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
#' \preformatted{
#' ## THIS FIELD SHOULD NOT BE MODIFIED BY THE USER OR OTHER CODE
#' ## Doing so will probably not cause problems, but I can't conceive
#' ## how it might help anything, either
#' }
#'
#' @return A data.table
#'
#' @seealso \link{reset}, \link{filterSummary}, \link{setFilters},
#'     \link{appliedFilters},
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

#' Set Filters
#'
#' Internal AnnotatedMatrix field holding a character vector
#' describing filters that have been set
#'
#' @name setFilters
#'
#' @details
#'
#' Each operation applied to the matrix to filter out values will add
#' a new string to this field. The string will be both (somewhat)
#' human-intelligible, as well as being machine-parsable. The
#' intention of the field is to allow reapplication of some or all of
#' the filters, as well as serialize the filters in an ASCII format
#' that can be stored with output of analyses.
#'
#' \preformatted{
#' ## THIS FIELD SHOULD NOT BE MODIFIED BY THE USER OR OTHER CODE
#' ## Doing so will probably not cause problems, but I can't conceive
#' ## how it might help anything, either
#' }
#'
#' @return A character vector
#'
#' @seealso \link{appliedFilters}, \link{filterSummary}
#'
#' @examples
#'
#' # The example matrix includes some automatic filters:
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$setFilters
NULL

#' Matrix Level Values
#'
#' Internal AnnotatedMatrix field holding factor levels (names)
#'
#' @name lvlVal
#'
#' @details
#'
#' AnnotatedMatrix objects can be defined as a "pseudo-factor" by
#' setting factor level names that will be associated with numeric
#' values of the matrix. These level names are held in the
#' \code{lvlVal} field. If the field is populated the matrix will be
#' presumed to be a factor.
#'
#' \preformatted{
#' ## THIS FIELD SHOULD NOT BE MODIFIED BY THE USER OR OTHER CODE
#' ## DOING SO COULD LEAD TO CONFUSION / CRASHES
#' }
#' 
#' @return A character vector of level values
#'
#' @seealso \link{is.factor}, \link{levels}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$lvlVal
NULL

#' Changed Row and Column Names
#'
#' Internal AnnotatedMatrix fields tracking row and column names that
#' needed to be changed when building the object
#'
#' @name rowcolchanges
#' @aliases rowChanges colChanges
#' @method rowChanges AnnotatedMatrix
#' @method colChanges AnnotatedMatrix
#' 
#'
#' @details
#'
#' Row and column names on the matrix must be unique. Some input
#' formats do not contain this restriction. If any dimension names
#' need to be altered to force them to be unique, those changes will
#' be captured here. The values are named character vectors, where the
#' value of the vector represents the original name, and the names
#' represent the uniquely-renamed values.
#'
#' \preformatted{
#' ## THIS FIELD SHOULD NOT BE MODIFIED BY THE USER OR OTHER CODE
#' ## Doing so will probably not cause problems, but I can't conceive
#' ## how it might help anything, either
#' }
#' 
#' @return A named character vector
#'
#' @seealso \link{rNames}, link{cNames}, link{uniqueNames}
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
#' myObject$matObj( help=TRUE )
#' 
#' obj <- myObject$matObj( raw=FALSE, transpose=FALSE )
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

#' Row and Column Names
#'
#' AnnotatedMatrix object method to get/set/reorder the rownames or
#' colnames of the filtered matrix
#'
#' @name rowcolnames
#' @aliases rNames cNames
#' @method rNames AnnotatedMatrix
#' @method cNames AnnotatedMatrix
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$rNames( help=TRUE )
#' 
#' rn <- myObject$rNames( raw=FALSE )
#' cn <- myObject$cNames( raw=FALSE )
#'
#' myObject$rNames( newRowNames, reason=NA )
#' myObject$cNames( newColNames, reason=NA )
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
#'     NULL. Optional human-readable reason for why the alteration was
#'     made, will be recorded in \link{filterLog} and used to
#'     structure \link{filterSummary}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character vector of row names
#'
#' @seealso \link{reset}, \link{rowChanges}, \link{colChanges}
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
#' myObject$reset( help=TRUE )
#' 
#' myObject$reset()
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
#' myObject$autoFilter( help=TRUE )
#'
#' myObject$autoFilter( recursive=TRUE, verbose=TRUE )
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
#'     \link{appliedFilters}, \link[ParamSetI]{showParameters}
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

#' Detail Zeroed Rows and Columns
#'
#' Internal AnnotatedMatrix object method to tally filtered rows and cols
#'
#' @name DOTdetailZeroedRowCol
#' @aliases .detailZeroedRowCol
#' @method .detailZeroedRowCol AnnotatedMatrix
#' 
#' @details
#' 
#' Internal method, should not be called directly
#' 
#' Takes a SparseMatrix and a (previously calculated) logical vector
#' of 'failed' i+j triples, and determines which (if any) rows and
#' columns went from populated to unpopulated (all zeroes). Adds
#' entries to \link{filterLog} to record them
NULL

#' Filter by Score
#'
#' AnnotatedMatrix object method to filter cells by their numeric value
#'
#' @name filterByScore
#' @method filterByScore AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$filterByScore( help=TRUE )
#'
#' myObject$filterByScore( min=NA, max=NA, filterEmpty=FALSE, reason=NA )
#' }
#' 
#' Updates the \link[=matrixUse]{used matrix} to remove (zero-out)
#' cells that fail a numeric test. Like all filters, will not alter
#' \link{matrixRaw}.
#'
#' @param min Optional minimum score. Cells that are below this value
#'     will be set to zero
#' @param max Optional maximum score. Cells that are above this value
#'     will be set to zero
#' @param filterEmpty Default FALSE. If true, then rows or columns
#'     that are empty after filtering will be removed from the matrix
#' @param reason Default NA. Optional human-readable reason for why
#'     the alteration was made, will be recorded in \link{filterLog}
#'     and used to structure \link{filterSummary}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return Invisibly, an integer vector of newly zeroed counts for
#'     (cells, rows, cols).
#'
#' @seealso \link{filterLog}, \link{filterSummary},
#'     \link{appliedFilters}
#'
#' @examples
#'
#' # The example matrix is factorized, but we can still fitler it by
#' # the numeric factor levels
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' 
#' # Check the levels:
#' s2e$levels()
#' 
#' # How many populated rows currently?
#' sum( s2e$populatedRows() )  # 151
#' 
#' # Let's only keep official symbols
#' s2e$filterByScore( min=4, reason="Focus on official symbols only" )
#' # We now have a more constrained matrix:
#' sum( s2e$populatedRows() )  # 28
#' s2e$filterSummary()
#' 
NULL

#' Filter by Factor Level
#'
#' AnnotatedMatrix object method to filter cells by factor level assignment
#'
#' @name filterByFactorLevel
#' @method filterByFactorLevel AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$filterByFactorLevel( help=TRUE )
#'
#' myObject$filterByFactorLevel(x, keep=TRUE, ignore.case=TRUE,
#'                              filterEmpty=FALSE)
#' }
#' 
#' Updates the \link[=matrixUse]{used matrix} to remove (zero-out)
#' cells that do not have the requested factor levels assinged. Like
#' all filters, will not alter \link{matrixRaw}.
#'
#' If the matrix is not factorized, an error will be reported, no
#' action will be taken, and NA will be invisibly returned.
#'
#' @param x Required, a list of factor levels. Can be either the
#'     factor level names (character) or integer level assignments.
#' @param keep Default TRUE, which will presume that \code{x}
#'     represents levels that should be kept. If FALSE, then levels in
#'     \code{x} will be removed.
#' @param ignore.case Default TRUE, which will match factor levels
#'     regardless of case.
#' @param filterEmpty Default FALSE. If true, then rows or columns
#'     that are empty after filtering will be removed from the matrix
#' @param reason Default NA. Optional human-readable reason for why
#'     the alteration was made, will be recorded in \link{filterLog}
#'     and used to structure \link{filterSummary}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return Invisibly, an integer vector of newly zeroed counts for
#'     (cells, rows, cols).
#'
#' @seealso \link{filterLog}, \link{filterSummary},
#'     \link{appliedFilters}
#'
#' @examples
#'
#' # The example matrix has some highly atypical symbols. To look at
#' # them, we will need to revert filters, since the file has
#' # autofilters that exclude them from the start:
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$reset()
#'
#' # Check the levels to make sure we recall what we want:
#' s2e$levels()
#'
#' s2e$filterByFactorLevel(c("Unknown"), keep=TRUE, filterEmpty=TRUE,
#'                         reason="Inspecting odd symbols")
#' s2e$filterSummary()
#'
#' # We set filterEmpty to TRUE in order to more easily inspect the
#' # matrix iteslf:
#' s2e$matrixUse
NULL

#' Filter by Metadata
#'
#' AnnotatedMatrix object method to filter rows or cells by metadata
#'
#' @name filterByMetadata
#' @method filterByMetadata AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$filterByMetadata( help=TRUE )
#'
#' myObject$filterByMetadata(key, val, MARGIN=NULL, type="like", op='or',
#'                           reason=NA)
#' }
#'
#' Updates the \link[=matrixUse]{used matrix} to remove (zero-out) all
#' rows or columns that have metadata failing to match your
#' request. Like all filters, will not alter \link{matrixRaw}.
#'
#' Note that this is not a cell-level filter, since metadata is
#' associated with the dimensions of the matrix, not the "connections"
#' between dimensions. It will thus remove entire rows or
#' columns. However, it will still tally the number of newly-zeroed
#' cells generated by the filter operation.
#'
#' @param key Required, the name of the key (tag / column) to match
#'     against. Will throw an error if not a single value
#' @param val Required, the value to match against. Can be multiple
#'     values.
#' @param MARGIN Default NULL, which will look for matches in both
#'     rows in columns. Alternatively can be 1, 2, "row" or "col".
#' @param keep Default TRUE, which will keep any matching entry,
#'     discarding (zeroing out) the others. If FALSE, then rows/cols
#'     with metadata which does not match will be removed.
#' @param type Default 'like', controls the kind of matching. 'like'
#'     will allow a fragmentary match of your value to metadata (both
#'     "left" and "right" wildcards), 'regexp' will treat your input
#'     as a regular expression, and "equal" will require an exact
#'     match. Additionaly, the text "case" can be added to any of the
#'     types to force a case-sensitive match; For example,
#'     "equal case" will require an exact, case-sensitive match.
#' @param op Default 'or', controls how filters resolve two or more
#'     \code{val} entries. 'or' will count a metadata entry as
#'     matching if any of your \code{val} entries match, while 'and'
#'     requires all to match. Additionally, you can use aliases 'any'
#'     (='or') or 'all (='and').
#' @param filterEmpty Default FALSE. If true, then rows or columns
#'     that are empty after filtering will be removed from the matrix
#' @param reason Default NA. Optional human-readable reason for why
#'     the alteration was made, will be recorded in \link{filterLog}
#'     and used to structure \link{filterSummary}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return Invisibly, an integer vector of newly zeroed counts for
#'     (cells, rows, cols).
#'
#' @seealso \link{filterLog}, \link{filterSummary},
#'     \link{appliedFilters}, \link{metadata}
#'
#' @examples
#' 
#' # The example matrix has some deprecated genes. These are already
#' # removed by the default auto filters, so we will reset the matrix
#' # before looking for them
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$reset()
#'
#' # This tag is used in the Description field to highlight deprecated loci:
#' depTag <- "{Deprecated}"
#' s2e$filterByMetadata("description", depTag, keep=TRUE, filterEmpty=TRUE,
#'                         reason="Inspecting deprecated loci")
#' s2e$filterSummary()
#'
#' # We set filterEmpty to TRUE in order to more easily inspect the
#' # matrix iteslf:
#' s2e$matrixUse
#'
#' # Verify that we found what we were looking for:
#' s2e$metadata( colnames(s2e$matrixUse), "Description" )
NULL


#' Number of Non-zero Cells
#'
#' AnnotatedMatrix object method to return a count of non-zero cells
#' in the used matrix
#'
#' @name nnZero
#' @method nnZero AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$nnZero( help=TRUE )
#'
#' myObject$nnZero( obj=NULL )
#' }
#'
#' For sparse matrices, the non-zero values represent the data, which
#' in AnnotatedMatrix are typically used to indicate connections of
#' some sort between the two dimensions of the matrix.
#'
#' @param obj Default NULL, which will recover the matrix using
#'     \link{matObj}. Optionally can be a user-provided sparse Matrix.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#' @param ... Passed to \link{matObj} (for example, \code{raw})
#'
#' @return An integer
#'
#' @seealso \link{counts}, \link{populated}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e
#' # The example matrix begins filtered - inspecting the object should show:
#' 
#' # 180 non-zero cells (1.92%)
#' # [Raw Matrix] : 167 x 67 (90.4%x92.5%), 208 non-zero (1.86%; 86.54% survive filters)
#'
#' # That is, in the filtered state 180 cells are non-zero. Confirm:
#' s2e$nnZero()           # 180
#' # We can see the total 'raw' connections by passing the raw argument:
#' s2e$nnZero( raw=TRUE ) # 208
#' # Or we can reset the filters and check again:
#' s2e$reset()
#' s2e$nnZero()           # 208
NULL

#' Row and Column Counts
#' 
#' AnnotatedMatrix object methods to count non-zero entries for rows or columns
#'
#' @name counts
#' @aliases rCounts cCounts
#' @method rCounts AnnotatedMatrix
#' @method cCounts AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$rCounts( help=TRUE )
#' myObject$cCounts( help=TRUE )
#'
#' myObject$rCounts( obj=NULL )
#' myObject$cCounts( obj=NULL )
#' }
#'
#' For \code{$rCounts()}, return the counts of non-zero columns for
#' each row. Same for \code{$cCounts()}, which will return the counts
#' of non-zero rows for each column.
#'
#' A value of zero indicates that the row (or column) does not have
#' any "connections" to the other dimension. A value of 1 indicates a
#' single unique connection.
#'
#' @param obj Default NULL, which will recover the matrix using
#'     \link{matObj}. Optionally can be a user-provided sparse Matrix.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#' @param ... Passed to \link{matObj} (for example, \code{raw})
#'
#' @return A named vector of integers, which will be as long as the
#'     number of rows (for rCounts) or columns (for cCounts). Names
#'     are taken from the dimension names.
#'
#' @seealso \link{populated}, \link{nnZero}
#'
#' @examples
#'
#' ## Example matrix has gene symbols in rows, gene IDs in columns
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' symCount <- s2e$rCounts()
#'
#' # The matrix is loaded with an automatic filter. Which symbols no
#' # longer connect to any IDs?
#' names( symCount[ symCount == 0 ] )
#'
#' # Which symbols are highly ambiguous (assigned to many genes)?
#' symCount[ symCount > 3 ]
#' # (Note that this matrix is case-sensitive, and includes "p40" and "P40")
#'
#' idCount <- s2e$cCounts()
#' # Which genes have a lot of symbols assigned to them?
#' idCount[ idCount > 7 ]
#' # What are they?
#' s2e$map( names(idCount[ idCount > 7 ]) )
NULL

#' Populated Rows and Columns
#' 
#' AnnotatedMatrix object methods to flag non-zero rows and columns
#'
#' @name populated
#' @aliases populatedRows populatedCols
#' @method populatedRows AnnotatedMatrix
#' @method populatedCols AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$populatedRows( help=TRUE )
#' myObject$populatedCols( help=TRUE )
#'
#' myObject$populatedRows( obj=NULL )
#' myObject$populatedCols( obj=NULL )
#' }
#'
#' These methods are just simple wrappers for the
#' \link[=counts]{counts methods}, which convert
#' their output to logical using \code{!= 0}. TRUE values indicate at
#' least one connection to the other dimension.
#'
#' @param obj Default NULL, which will recover the matrix using
#'     \link{matObj}. Optionally can be a user-provided sparse Matrix.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#' @param ... Passed to \link{matObj} (for example, \code{raw})
#'
#' @return A named logical vector, which will be as long as the
#'     number of rows (for rCounts) or columns (for cCounts). Names
#'     are taken from the dimension names.
#'
#' @seealso \link{counts}, \link{nnZero}
#'
#' @examples
#'
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' # Reset filters, then toss in a metadata filter to make a smaller matrix
#' s2e$reset()
#' s2e$filterByMetadata("Description","superfamily")
#'
#' # What genes and loci still have connections after this?
#' popRows <- s2e$populatedRows()
#' popRows[ popRows ]
#'
#' popCols <- s2e$populatedCols()
#' s2e$metadata(names(popCols[ popCols ]), c("Symbol", "Description"))
#' 
NULL

#' Remove Empty Rows and Columns
#' 
#' AnnotatedMatrix object methods to remove rows or columns that are all zeroes
#'
#' @name removeEmpty
#' @aliases removeEmptyRows removeEmptyCols
#' @method removeEmptyRows AnnotatedMatrix
#' @method removeEmptyCols AnnotatedMatrix
#' @method removeEmpty AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$removeEmptyRows( help=TRUE )
#' myObject$removeEmptyCols( help=TRUE )
#' myObject$removeEmpty( help=TRUE )
#'
#' myObject$removeEmptyRows( reason=NA )
#' myObject$removeEmptyCols( reason=NA )
#' myObject$removeEmpty( reason=NA )
#' }
#'
#' AnnotatedMatrix can work with matrices that have entirely "empty"
#' rows or columns; For sparse matrices these would correspond to
#' those that are all zeroes.
#'
#' There are advantages to leaving such entries in the matrix: In
#' \code{$map()}, the presence of an empty row/col can be used to
#' differentiate between "This input no longer has connections" versus
#' "This input is unrecognized". Additionally, empty entries allow for
#' name-based access without throwing errors.
#'
#' However, in some cases it may be desirable to entirely remvove rows
#' or columns. In particular, crossproducts between matrices requires
#' the length of the crossed dimension to match.
#'
#' These operations, like all filters, will only affect the
#' \link[=matrixUse]{used matrix}. They will not generate new
#' \link{filterLog} entries, however. Rather, any prior fitlers that
#' resulted in a row or column becoming fully zeroed-out will have
#' noted the event. The action will create an entry in
#' \link{setFilters}, however.
#'
#' \code{$removeEmpty()} simply wraps both \code{$removeEmptyRows()}
#' and \code{$removeEmptyCols()}
#'
#' @param reason Default NA, optional human-readable text that will be
#'     included with the \link{setFilters} entry.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return Invisibly, a character vector of removed names
#' 
#' @seealso \link{counts}, \link{nnZero}, \link{populated}
#'
#' @examples
#' 
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' # The example matrix has some filters already applied:
#' rc <- s2e$rCounts()
#' length(rc)    # Total number of rows = 167
#' rc[ rc == 0 ] # Rows that no longer connect to columns
#'
#' gone <- s2e$removeEmptyRows()
#' gone                 # Removed names
#' rce <- s2e$rCounts()
#' length(rce)          # -> 151
#' rce[ rce == 0 ]      # No more empty rows
NULL


#' Map
#' 
#' AnnotatedMatrix object method to convert names on one dimension to the other
#' 
#' @name map
#' @method map AnnotatedMatrix
#' 
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$map( help=TRUE )
#'
#' myObject$map(input=NULL, via=NULL, format="data.frame", ignore.case=TRUE,
#'              keep.best=FALSE, column.func=max, collapse=NULL,
#'              collapse.name=NULL, collapse.token=',',
#'              collapse.score=NULL, collapse.factor=NULL,
#'              integer.factor=FALSE, add.metadata=TRUE, warn=TRUE,
#'              append.to=NULL, append.col=1L, help=FALSE )
#' }
#'
#' This is the primary method for AnnotatedMatrix. It takes the input,
#' matches it to one of the dimensions (chosen either automatically or
#' specifically) and then reports "connections" to the other
#' dimension. The results can be returned in a variety of formats.
#'
#' The method is alert for any non-unique mappings, and will record
#' such cases in attributes of the returned value as well as by alert
#' to the terminal. The following non-unique categories are noted:
#'
#' \preformatted{
#'  Dup.In   : IDs that were present twice or more in input
#'             (possibly after case removal)
#'  Dup.Mat  : IDs that were present twice or more in matrix
#'             (always after case removal, since matrix names are unique)
#'  Mult.In  : Input IDs that generated multiple output values
#'  Mult.Out : Output IDs that were generated from multiple inputs
#'  Unmapped : Input ID is also in the matrix, but does not have a target
#'             with non-zero score. Score will be zero
#'  Unknown  : Input ID could not be matched to any in the matrix. Score
#'             will be \code{NA}
#' }
#'
#' These categories also correspond to attributes of the same names
#' that are attached to the returned value. In addition, a "Notes"
#' attributes contains the above descriptions to aid in remembering
#' the distinctions between categories.
#'
#' Note also that there are two "not mapped" scores. \code{NA}
#' indicates that the requested input value was not found in the
#' matrix, either because it was never there at all, or because the
#' row/col that held it was removed (through \link{removeEmpty} or
#' \link{rNames} / \link{cNames}). Alternatively, a value of \code{0}
#' indicates that the name is recognized, but now occupies a fully
#' zeroed-out row or column (presumably through applied
#' filters). These two distinct values can be helpful in
#' troubleshooting why some input values are not producing output.
#'
#' @param input Default NULL. The method requires input to
#'     function. It can be provided here as a character vector of IDs,
#'     or if \code{append.to} has provided a data.frame, then it will
#'     be taken from the column in that object specified by
#'     \code{append.col}
#' @param via The matrix dimension that \code{input} is from. The
#'     default is NULL, which will compare your input to both rows and
#'     columns, and chose the dimension with the most matches
#'     (defaulting to 'row' in the event of equal matches)
#' @param format Default 'data.frame', specifies the output
#'     format. Can also be 'vector', which will return a named
#'     character vector, with values corresponding to \code{$Output}
#'     and names as \code{$Input}
#' @param ignore.case Default TRUE, which ignores the capitilazation
#'     of IDs in both columns and rows
#' @param keep.best Default FALSE. If TRUE, then only the top-scored
#'     cell(s) will be kept for each input conversion
#' @param column.func Default \code{max}.  If ignore.case is true, it
#'     is possible that an input ID can match multiple matrix IDs. In
#'     this case, multiple matching rows will be returned for one
#'     ID. \code{column.func} is applied to reduce this to a single
#'     row. The function should accept a numeric vector of scores as
#'     the first argument.
#' @param collapse Default NULL, which will cause every pairwise
#'     connection to be reported. If 'in', then the \code{$Input}
#'     column of the data.frame will be unique - any input value that
#'     results in multiple output values will result in the
#'     \code{$Output} IDs and \code{$Score} being 'collapsed' to a
#'     single value (see the \code{collapse.*} options below). A value
#'     of 'out' will do the same, but causes the \code{$Output} column
#'     to be unique, and \code{$Input} and \code{$Score} are
#'     collapsed.
#' @param collapse.name Default NULL, which will cause multiple names
#'     to be concatenated into a single value using
#'     paste(). Alternatively, a user function can be provided. This
#'     package also includes the crude utility function
#'     \link{takeLowestThing} which can be used here. The function
#'     should accept a character vector of names as the first
#'     argument.
#' @param collapse.token Default ',', a string used to concatenate
#'     collapsed IDs when using paste via \code{collapse.name=NULL}.
#'     An error will be thrown if it is not a single scalar string.
#' @param collapse.score Default NULL. Optional function to apply to
#'     the \code{$Score} column when two or more rows are being
#'     collapsed into one. The function should take a numeric vector
#'     as input and return a single numeric value. If NULL, the
#'     function will be mean(), unless the matrix is a pseudo-factor,
#'     in which case the object method \link{.autoLevel} will be
#'     used. This will generate new 'hybrid' factors as
#'     needed. \code{collapse.token} will be used as the string to
#'     join factors into a new level name.
#' @param integer.factor Default NULL, which will include both the
#'     \code{$Score} column (integer factor value) as well as a
#'     \code{$Factor} column (with level values as characters) in the
#'     output to be present instead. If TRUE then ONLY a \code{$Score}
#'     column (representing integer values, perhaps including new and
#'     likely-meaningless hybrid values from \link{.autoLevel}) will
#'     be added. If FALSE, then only the \code{$Factor} column will be
#'     present.
#' @param add.metadata Default TRUE, which will add all metadata
#'     columns that have at least some information associated with the
#'     \code{$Ouput} column. FALSE will prevent adding metadata, and a
#'     character vector will add those specific columns (it is up to
#'     the user to confirm that requested columns exist in the
#'     metadata store)
#' @param warn Default TRUE, which will show warning text if matches
#'     failed to be made for the input. This information is also
#'     always captured in attributes attached to the returned
#'     data.frame
#' @param append.to Default NULL. If a data.frame-compliant object, it
#'     will become the return value, with the columns generated by
#'     this function being added on. If \code{collapse='out'} then the
#'     \code{$Output} column will be used to merge, otherwise the
#'     \code{$Input} column will be utilized. The merge column from
#'     the provided data.frame is set with \code{append.col}
#' @param append.col The column to use in \code{append.to} for
#'     merging. Default is \code{1L}, can provide another column
#'     number or name. If \code{input} was not specified, then the
#'     contents of this column will be taken as input.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#' 
#' @return If \code{format} is 'data.frame', a data.frame. If
#'     \code{format} is 'vector', an SVG map of the globe, with output
#'     values represented as minor bodies of water or occasionally
#'     rivers. Just kidding, you'd get a (named) vector.
#' 
#' @seealso \link{.autoLevel}, \link{metadata}
#'
#' @examples
#' 
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' 
#' # I have some gene symbols, what Entrez Gene IDs do they map to?
#' syms <- c("FOXD3", "NBP", "AR")
#' s2e$map( syms )
#'
#' # More messy than I wanted, "NBP" and "AR" are not uniquely
#' # mapping. Let's filter the matrix to only include official symbols
#' s2e$filterByFactorLevel("Official", keep=TRUE)
#' s2e$map( syms )
#'
#' # Now we've lost mappings to NBP. Let's try again, but now use keep.best
#' s2e$reset() # Clear the filters (will also clear an automatic filter)
#' s2e$map( syms, keep.best=TRUE )
#'
#' # The two NBP mappings are both unofficial. We could work with the
#' # data as-is, or we could collapse on input to assure that each row
#' # is unique:
#' s2e$map( syms, keep.best=TRUE, collapse='in' )
#'
#' # We can also merge these data into an existing structure. Say we
#' # have some gene-linked research data:
#' evilCRO <- data.frame(
#'      Client=c("Evil Inc.", "Quite Evil", "Totally Evil"),
#'      Mutation=c("FOXD3", "NBP", "AR"),
#'      Phenotype=c("Eye lasers", "Adamantine claws", "Bunny whiskers"),
#'      Contained=c(TRUE,TRUE,FALSE),
#'      stringsAsFactors=FALSE)
#'
#' # We can generate an integrated data.frame with append.to. Don't
#' # forget to specify the merge column with append.col if it's not
#' # the first column.
#' evilLoci <- s2e$map(append.to=evilCRO, append.col="Mutation",
#'                     keep.best=TRUE, collapse='in' )
#' evilLoci
#'
#' # For large results the non-unique events may not be
#' # transparent. Check the attributes to identify non-unique IDs
#' attr(evilLoci, 'Mult.In')
#' 
NULL

#' Auto Level
#'
#' Internal AnnotatedMatrix object method to generate new hybrid factor levels
#'
#' @name DOTautoLevel
#' @aliases .autoLevel
#' @method .autoLevel AnnotatedMatrix
#' 
#' @details
#' 
#' Internal method, should not be called directly
#' 
#' Takes a set of factor levels (character names), combines them via
#' paste, generates a new factor level if needed, and returns the
#' factor index. Will extend the \code{lvlVal} internal field as
#' needed.
#'
#' This function is used by \link{map} to generate synthetic factor
#' levels when \code{collapse} causes multiple rows to be merged into
#' a single one.
#'
#' @param vals Required, vector of one or more factor levels
#' @param sep Default ',', the token used to separate multiple level
#'     names
#' @param decreasing Default TRUE, which will sort factor values from
#'     largest to smallest (in order to put the 'best'-scored factor
#'     at the front of the list)
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return An integer representing the factor level
#'
#' @seealso \link{map}
NULL

#' Filter Details
#'
#' Internal AnnotatedMatrix object method to add filtered names to the log
#'
#' @name DOTfilterDetails
#' @aliases .filterDetails
#' @method .filterDetails AnnotatedMatrix
#' 
#' @details
#' 
#' Internal method, should not be called directly
#' 
#' Helper method to add new filter entries to \link{filterLog}.
#'
#' @param id Required, a vector of one or more IDs (row or column
#'     names) that have been filtered out.
#' @param type Optional type of filter, generally Row, Col or Cell
#' @param metric Optional human-readable criteria for the filter, such
#'     as "Score > 3.14"
#' @param reason Optional human-readable rationale for why the filter
#'     was applied, such as "Remove low-confidence assignments"
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The \link{filterLog} filed (data.table) itself
#'
#' @seealso \link{filterLog}, \link{filterSummary},
#'     \link{appliedFilters}
NULL

#' Add Applied Filter
#'
#' Internal AnnotatedMatrix object method to extend the list of applied filters
#'
#' @name DOTaddAppliedFilter
#' @aliases .addAppliedFilter
#' @method .addAppliedFilter AnnotatedMatrix
#' 
#' @details
#' 
#' Internal method, should not be called directly
#' 
#' Helper method to add new filter entries to \link{setFilters}.
#'
#' @param key Required, a key indicating the filter type, eg 'SCORE' or ''
#' @param val Required, text defining the filter parameters
#' @param com Default NULL, optional comment (reason) for the filter
#' @param pre Default NULL, optional text (eg operator) to go before the values
#' @param pro Default NULL, text (eg modifiers) to go after the values
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A single character string representing the fitler text
#'     added to \link{setFilters}
#'
#' @seealso \link{filterLog}, \link{filterSummary},
#'     \link{appliedFilters}
NULL

#' Applied Filters
#'
#' AnnotatedMatrix object method to show/add applied filters
#'
#' @name appliedFilters
#' @method appliedFilters AnnotatedMatrix
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$appliedFilters( help=TRUE )
#'
#' myObject$appliedFilters( new=NULL )
#' }
#'
#' Each time a method is run that filters the matrix, the criteria
#' used are serialized as a text snippet, and the snippet is added to
#' the \link{setFilters} vector. These snippets can be used as a
#' record of filtering operations. In addition, some or all of them
#' can be passed to this function to re-run the filters
#'
#' (LOADING = ASPIRATIONAL FEATURE, NOT YET AVAILABLE).
#'
#' @param new Default NULL, optional character vector of new filters
#'     to apply. NOT YET IMPLEMENTED
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The \link{setFilters} field, a character vector
#'     added to \link{setFilters}
#'
#' @seealso \link{setFilters}, \link{filterByScore},
#'     \link{filterByFactorLevel}, \link{filterByMetadata},
#'     \link{rNames}, \link{cNames}, \link{reset},
#'     \link{filterSummary},
#'
#' @examples
#' 
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' ## The sample file has automated filters applied to it:
#' s2e$appliedFilters()
NULL

#' Filter Summary
#'
#' AnnotatedMatrix object method to summarize filters applied to the matrix
#'
#' @name filterSummary
#' @method filterSummary AnnotatedMatrix
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$filterSummary( help=TRUE )
#'
#' myObject$filterSummary( reason=TRUE )
#' }
#'
#' Returns a data.frame that has been classed to allow a pretty-print
#' summary of the filters applied and the number of rows or columns
#' that have been 'removed' by each filter.
#'
#' @param reason Default TRUE, which will cluster the filter summary
#'     by reason and type. If FALSE, then only type will be used to
#'     cluster the summary.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A data.frame classed as 'mapFilter'. Will be pretty-printed
#'     with \link{print.mapFilter}.
#'
#' @seealso \link{setFilters}, \link{filterByScore},
#'     \link{filterByFactorLevel}, \link{filterByMetadata},
#'     \link{rNames}, \link{cNames}, \link{reset},
#'     \link{appliedFilters},
#'
#' @examples
#' 
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' fs <- s2e$filterSummary( )
#' fs
#' str(fs)
NULL

#' Matrix is Factor
#'
#' AnnotatedMatrix object method indicating if the matrix is a factor or not
#'
#' @name is.factor
#' @method is.factor AnnotatedMatrix
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$is.factor( help=TRUE )
#'
#' myObject$is.factor( )
#' }
#'
#' \link{Matrix} matrices do not support factorized content, at least
#' at the time this package was developed. AnnotatedMatrix maintains
#' some metadata to allow integer matrices to be treated as factors
#' with named levels. This method will indicate if the \link{lvlVal}
#' field has ben set.
#'
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A single logical value, TRUE if \link{lvlVal} is defined
#'     with one or more levels, otherwise FALSE.
#'
#' @seealso \link{lvlVal}, \link{levels}
#'
#' @examples
#' 
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$is.factor( )
#' s2e$levels()
NULL

#' Matrix Levels
#'
#' AnnotatedMatrix object method returning the factor levels (names)
#'
#' @name levels
#' @method levels AnnotatedMatrix
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$levels( help=TRUE )
#'
#' myObject$levels( )
#' }
#'
#' If the matrix is flagged as being a factor (by having \link{lvlVal}
#' set), this function will return the factor levels (names).
#'
#' @param asFactor Default FALSE, in which case factor level names
#'     only will be returned. If TRUE, then the returned value will be
#'     a formal factor object.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return NULL if the matrix is not a factor, or a character vector
#'     of level names, or a factor if \code{asFactor} is TRUE.
#'
#' @seealso \link{lvlVal}, \link{is.factor}
#'
#' @examples
#' 
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$levels()
#' str( s2e$levels( asFactor=TRUE ) )
NULL

#' As GMT
#'
#' AnnotatedMatrix object method to represent matrix as GMT format
#'
#' @name as.gmt
#' @method as.gmt AnnotatedMatrix
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$as.gmt( help=TRUE )
#'
#' myObject$as.gmt( obj=NULL, transpose=FALSE, file=NULL, ...)
#' }
#'
#' Will represent the matrix as GMT format, a simple flat file format
#' that allows for basic metadata (names and descriptions).
#'
#' @param obj Default NULL, the matrix object to serialize. If NULL,
#'     will call \link{matObj} to recover the current filtered matrix.
#' @param transpose Default FALSE, which will treat rows as sets and
#'     columns as set members. TRUE will instead treat the columns as
#'     sets.
#' @param file Default NULL. If not NULL, will be treated as a file
#'     path and passed to \link{cat}.
#' @param ... Will be passed to \link{matObj}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The GMT text, invisibly if file is not NULL
#'
#' @seealso \href{http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{GMT Format at Broad Institute}
#'
#' @examples
#' 
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' ## Generate GMT data for the symbols associated for each gene in
#' ## the example:
#' s2e$as.gmt( transpose=TRUE )
NULL

#' Read Matrix Data
#'
#' Internal AnnotatedMatrix object method to parse file contents
#'
#' @name DOTreadMatrix
#' @aliases .readMatrix
#' @method .readMatrix AnnotatedMatrix
#' 
#' @details
#' 
#' Internal method, should not be called directly
#' 
#' This method will take the \link{file} field as a file path and
#' attempt to parse it as a matrix. Once parsed, the method will
#' create a serailized \link[=readRDS]{RDS} object. Future attempts to
#' parse the matrix file will recognize the RDS object and quickly
#' parse it instead.
#'
#' @param format Default "" (empty string). Optional format name, if
#'     empty the format will be deduced by the file suffix.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A list containing the \link{Matrix} object, metadata as a
#'     \link{data.table}, \link{rowChanges} and \link{colChanges},
#'     levels as \link{lvlVal}, and any default parameters parsed from
#'     the file.
#'
#' @seealso The individual format parsing functions:
#'     \link{parse_MatrixMarket_file}, link{filesToMTX},
#'     link{parse_ListOfLists_file}, link{parse_Text_file},
#'     link{parse_GMT_file}
NULL

#' Metadata Keys
#'
#' AnnotatedMatrix object method to return all metadata keys
#'
#' @name metadata_keys
#' @method metadata_keys AnnotatedMatrix
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$metadata_keys( help=TRUE )
#'
#' myObject$metadata_keys( )
#' }
#'
#' Return all keys (tag names) from the internal metadata table. The
#' \code{id} column will be excluded from the result.
#'
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character vector of metadata column names
#'
#' @seealso \link{metadata}, \link{matrixMD}
#'
#' @examples
#' 
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' ## Two keys assigned to the example matrix:
#' s2e$metadata_keys()
NULL

#' Metadata Keys
#'
#' AnnotatedMatrix object method to recover metadata assigned to rows or columns
#'
#' @name metadata
#' @method metadata AnnotatedMatrix
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$metadata( help=TRUE )
#'
#' myObject$metadata(id=NULL, key=NULL, na.rm=TRUE, drop=TRUE, verbose=TRUE )
#' }
#'
#' Query metadata associated with rows and/or columns, based on IDs
#' and/or keys
#'
#' @param id Default NULL. Optional vector of specific IDs to query
#' @param key Default NULL. Optional vector of keys / tags to query
#' @param na.rm Default TRUE, which will cause rows or columns that
#'     are entirely NA to be removed.
#' @param drop Default TRUE, which will cause a named vector to be
#'     returned if only one column of metadata is present. Otherwise a
#'     data.table will be returned.
#' @param verbose Default TRUE, which will warn about certain issues
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#' 
#' @return If drop=TRUE and zero or one columns of metadata a present,
#'     a named character vector. Otherwise a data.table
#'
#' @seealso \link{metadata_keys}, \link{matrixMD}
#'
#' @examples
#' 
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$metadata(c("LOC9588","LOC7038","LOC3921"))
#' s2e$metadata(key="Symbol")
NULL

#' Field Descriptions
#'
#' Internal AnnotatedMatrix object method providing descriptive text
#' to object fields
#'
#' @name DOTfieldDescriptions
#' @aliases .fieldDescriptions
#' @method .fieldDescriptions AnnotatedMatrix
#' 
#' @details
#' 
#' Internal method, should not be called directly
#' 
#' Helper method that adds descriptive attributes to the fields held
#' by the object. This is an attempt to make the object
#' self-documenting. Two fields are added: "Description", which
#' provides a one-sentance description of the field, and "Help", which
#' is a copy-and-pastable \link{help} command to recover the man page
#' for the field
#'
#' @param update Default TRUE, which will modify the internal object
#'     field to have the two attached objects.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A list, with field variable names as names, and descriptive
#'     text as values.
#' @seealso \link{help}
NULL

