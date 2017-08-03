sfSep <- ' || ' # Token for separating text while recording filters

#' Annotated Matrix
#'
#' Annotated sparse matrix for capturing query lists, identifier
#' mappings and ontology lookups
#'
#' @details
#'
#' \preformatted{
#' ## Object creation
#' AnnotatedMatrix( help=FALSE ) # Show this help
#'
#' myAnnMat <- AnnotatedMatrix( file=NA, params=NA, autofilter=TRUE, ...)
#'
#' myAnnMat$help()                 # High-level help
#' myAnnMat$ANYMETHOD( help=TRUE ) # Detailed help for all methods
#' }
#'
#' AnnotatedMatrix is a heavy wrapper around the Matrix package, in
#' particular the sparse matrix dgTMatrix class. It was built with an
#' eye toward using the MatrixMarket file format, an efficient but
#' unfortunately simple (no structured annotation) file format. This
#' package will parse information from the comments ('%') of
#' MatrixMarket files to generate rich objects for enhanced filtering
#' and reporting.
#'
#' The driving motivation for the package is to support enrichment
#' analysis (Fisher's exact test) in a reproducible research
#' framework. The package is also very well suited for identifier
#' mapping (for example, gene symbols to gene accessions).
#'
#' @field file Path to file the matrix was loaded from
#' @field fromRDS Logical, true if the loaded file was an RDS object
#' @field log EventLogger object holding log (activity) entries
#' @field matrixRaw The Matrix object as loaded from the file, will
#'     not be altered by the object.
#' @field matrixUse The Matrix object after manipulation by any
#'     applied filters
#' @field matrixMD data.table holding metadata associated with the
#'     matrix
#' @field lvlVal Character array of level names for factor matrices
#' @field filterLog data.frame storing filtering events that transpire
#'     during the pruning of matrices prior to analysis.
#' @field setFilters Human- and machine-readable character vector of
#'     filters that have been applied
#' @field rowChanges Named character vector of any row names that
#'     needed changing. Values are the original name, names are the
#'     names after processing with make.names() (if valid = TRUE) or
#'     make.unique() if (valid = FALSE)
#' @field colChanges As per rowChanges, but for column names
#'
#' @importFrom methods new setRefClass
#' @importFrom utils read.table
#' @importFrom data.table data.table as.data.table setkeyv
#' @importFrom dplyr count
#' @import Matrix
#' @importFrom CatMisc is.def is.something is.empty.field methodHelp
#' @import ParamSetI
#' @importClassesFrom EventLogger EventLogger
#'
#' @examples
#'
#' ## In most circumstances you will provide only a file path, but
#' ## optionally a list of parameters and a flag to turn off autoFilter
#' ## can be provided as well.
#' 
#' ## A toy (small) symbol-to-gene mapping matrix is provided with the
#' ## package:
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' 
#' s2e$help()         # Compact reminder of fields and methods
#' s2e                # Summary of loaded object
#' s2e$map(help=TRUE) # Help for the map() method
#' AnnotatedMatrix(help=TRUE) # display this help
#' 
#' ## Load a toy symbol-to-gene mapping matrix and use it to convert
#' ## some genes to Entrez Gene IDs
#' demo("geneSymbolMapping", package="AnnotatedMatrix", ask=FALSE)
#'
#' ## Work with different file formats
#' demo("fileFormats", package="AnnotatedMatrix", ask=FALSE)
#'
#' ## Different approaches to handling non-unique mappings:
#' demo("managingMultiplicity", package="AnnotatedMatrix", ask=FALSE)
#'
#' ## Working with metadata
#' demo("workingWithMetadata", package="AnnotatedMatrix", ask=FALSE)
#' 
#' @export AnnotatedMatrix
#' @exportClass AnnotatedMatrix
#' 

AnnotatedMatrix <-
    setRefClass("AnnotatedMatrix",
                fields = list(
                    file       = "character",
                    fromRDS    = "logical",
                    matrixRaw  = "dgTMatrix",
                    matrixUse  = "ANY", # dgTMatrix
                    matrixMD   = "data.table",
                    filterLog  = "data.table",
                    setFilters = "character",
                    lvlVal     = "character",

                    ## If row or column names need to be remapped:
                    rowChanges = "character",
                    colChanges = "character"
                    ),
                contains = c("EventLogger", "ParamSetI")
                )


AnnotatedMatrix$methods(
    
    initialize = function(file=NA, params=NA, autofilter=TRUE,
                          help=FALSE, ... ) {
        "Create a new AnnotatedMatrix object; Invoke with AnnotatedMatrix(...)"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }
        callSuper(...)
        if (!CatMisc::is.def(file))
            err("AnnotatedMatrix must define 'file' when created", fatal = TRUE)
        file    <<- file
        fromRDS <<- FALSE
        defineParameters("
Name        [character] Short Name assigned to the matrix
Description [character] Description for the matrix
ScoreDesc   [Character] Describes what the matrix values (scores) represent
Source      [character] Primary source, presumably a URL
Authority   [character] The name of the authority responsible for the data

RowDim      [character] Name for the row dimension
ColDim      [character] Name for the column dimension
RowUrl      [character] Base URL for row names (%s placeholder for name)
ColUrl      [character] Base URL for column names (%s placeholder for name)

MinScore    [numeric] Minimum score recognized by $autoFilter()
MaxScore    [numeric] Maximum score recognized by $autoFilter()
KeepLevel   [character] List of preserved factor levels recognized by $autoFilter()
TossLevel   [character] List of discarded factor levels recognized by $autoFilter()
TossMeta    [character] Metadata value filter recognized by $autoFilter()

")
        setParamList( params=params, ... )
        .readMatrix( ... )
        reset()
        if (autofilter) autoFilter()
    },
    
    matObj = function( raw=FALSE, transpose=FALSE, help=FALSE) {
        "Return the current filtered state of the dgTMatrix sparse matrix"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        rv <- if (raw || !CatMisc::is.def(matrixUse)) {
                  matrixRaw } else { matrixUse }
        if (transpose) {
            ## Need to use Matrix's ::t function, AND need to manage
            ## dimension names - not sure if I am missing something
            ## here, but they did not copy over
            dims <- dimnames(rv)
            dnn  <- names(dims) # dimension names
            ## Awkward problems when presuming that both dimensions
            ## are named. Being cautious in setting up the new dimname list:
            smid <- list()
            if (is.def(dnn[2])) {
                smid[[ dnn[2] ]] <- dims[[2]]
            } else {
                smid[[ 1 ]] <- dims[[2]]
            }
            if (is.def(dnn[1])) {
                smid[[ dnn[1] ]] <- dims[[1]]
            } else {
                smid[[ 2 ]] <- dims[[1]]
            }
            rv <- Matrix::t( rv )
            dimnames(rv) <- smid
        }
        rv
    },

    rNames = function(new=NULL, raw=FALSE, reason=NA, help=FALSE) {
        "Get/set rownames for the matrix"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        obj   <- matObj(raw)
        if (is.null(new)) {
            ## Just asking for the current names
            return (rownames(obj))
        }
        ## User is also reodering names, or removing/adding some
        if (raw) err("$rNames() can not set new names while specifying raw=TRUE", fatal=TRUE)
        obj   <- matObj(raw=raw)
        now   <- rownames(obj)
        lost  <- setdiff(now, new)
        if (length(lost > 0)) {
            ## Some rows have been removed. But did they have any
            ## nonzero cells?
            counts <- rCounts( obj )[ lost ]
            counts <- counts[ counts > 0 ]
            if (length(counts) > 0) {
                ## Yes, at least one data-bearing row has been
                ## eliminated by this action
                .filterDetails(id=names(counts), type="Row", reason=reason,
                              metric="Row names reordered/reset")
            }
        }
        gain  <- setdiff(new, now)
        if (length(gain) > 0) {
            ## There are some new rownames, we will add as "empty"
            ## rows in the expected place (row order)

            ## The Matrix package has a problem with NA indices, so
            ## this does NOT work:
            ## obj <- obj[ match(new, now), , drop=FALSE ]

            ## Instead, extract the known rows:
            both <- intersect(new, now)
            obj  <- obj[ both, , drop=FALSE ]
            ## ... build a new empty matrix with the new ones.
            m <- Matrix(0, length(gain), ncol(obj), dimnames = list(gain, NULL))
            ## ... bind them into the known rows, and clean up the order:
            obj  <- rbind(obj, m)[ new, ]
            ## Both Matrix() and rbind() generate class 'dgCMatrix'
            ## objects. Coerce it back to dgTMatrix
            matrixUse <<- as(obj, "dgTMatrix")
            .filterDetails(id=gain, type="Row", reason=reason,
                          metric="Empty rows added")
        } else {
            ## No new rows, just a reshuffle and/or removal of existing names
            matrixUse <<- obj[ new, , drop=FALSE ]
        }
        .addAppliedFilter("ROWORDER", new, reason)
        ## Just return the provided value back to the user
        new
    },

    cNames = function(new=NULL, raw=FALSE, reason=NA, help=FALSE) {
        "Get/set colnames for the matrix"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        obj   <- matObj(raw)
        if (is.null(new)) {
            ## Just asking for the current names
            return (colnames(obj))
        }
        ## User is also reodering names, or removing/adding some
        if (raw) err("$cNames() can not set new names while specifying raw=TRUE", fatal=TRUE)
        now   <- colnames(obj)
        lost  <- setdiff(now, new)
        if (length(lost > 0)) {
            ## Some columns have been removed. But did they have any
            ## nonzero cells?
            counts <- rCounts( obj )[ lost ]
            counts <- counts[ counts > 0 ]
            if (length(counts) > 0) {
                ## Yes, at least one data-bearing column has been
                ## eliminated by this action
                .filterDetails(id=names(counts), type="Col", reason=reason,
                              metric="Col names reordered/reset")
            }
        }
        gain  <- setdiff(new, now)
        if (length(gain) > 0) {
            ## There are some new colnames, we will add as "empty"
            ## cols in the expected place (column order)

            ## The Matrix package has a problem with NA indices, so
            ## this does NOT work:
            ## obj <- obj[ , match(new, now), drop=FALSE ]

            ## Instead, extract the known columns:
            both <- intersect(new, now)
            obj  <- obj[ , both, drop=FALSE ]
            m <- Matrix(0, nrow(obj), length(gain), dimnames = list(NULL, gain))
            ## ... bind them into the known rows, and clean up the order:
            obj  <- cbind(obj, m)[ , new ]
            ## Both Matrix() and rbind() generate class 'dgCMatrix'
            ## objects. Coerce it back to dgTMatrix
            matrixUse <<- as(obj, "dgTMatrix")
            .filterDetails(id=gain, type="Col", reason=reason,
                           metric="Empty columns added")
        } else {
            ## No new columns, just a reshuffle and/or removal of existing names
            matrixUse <<- obj[ , new, drop=FALSE ]
        }
        .addAppliedFilter("COLORDER", new, reason)
        ## Just return the provided value back to the user
        new
    },

    reset = function( help=FALSE ) {
        "Reset the matrix to the 'raw' (unfiltered) state"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        matrixUse  <<- NULL
        setFilters <<- character()
        filterLog  <<- data.table(id   = character(), metric = character(),
                                  type = character(), reason = character(),
                                  key = "id")
        invisible(NA)
    },

    autoFilter = function( recursive=TRUE, verbose=TRUE, help=FALSE) {
        "Automatically apply filters defined in the parameters"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        "\\preformatted{
Will apply filters based on certain parameters. The method is called
automatically on object instantiation, but can be called manually as
well. Parameters will generally be set as part of the source data file.

  recursive - Default TRUE, which will result in autoFilter being
              recursively called until no further changes are made.
              This can be relevant with row (or column) count filters,
              where rows that pass on one iteration fail on another,
              due to zero-ed out cells that previously enabled a count
              to pass.
    verbose - Default TRUE, will print a message of how many cells have
              been affected.

Recognized parameters:
   MinScore   MaxScore
   KeepLevel  TossLevel
}"
        x <- 0
        min <- param("MinScore")
        if (CatMisc::is.something(min)) {
            x <- x + filterByScore(min=min, reason=attr(min,"comment"))
        }
        max <- param("MaxScore")
        if (CatMisc::is.something(max)) {
            x <- x + filterByScore(max=max, reason=attr(max,"comment"))
        }
        if (is.factor()) {
            kL <- param("KeepLevel")
            if (CatMisc::is.something(kL)) {
                x <- x + filterByFactorLevel(kL, keep=TRUE,
                                             reason=attr(kL,"comment"))
            }
            tL <- param("TossLevel")
            if (CatMisc::is.something(tL)) {
                x <- x + filterByFactorLevel(tL, keep=FALSE,
                                             reason=attr(tL,"comment"))
            }
        }
        meta <- param("TossMeta")
        if (CatMisc::is.something(meta)) {
            fbms <- CatMisc::parenRegExp('^(\\S+)\\s+(\\S+)\\s+(.+)', meta,
                                         unlist=FALSE)
            for (fbm in fbms) {
                v <- fbm[3]
                ## ToDo - smuggle in Row/Col at front of val
                x <- x + filterByMetadata(key=fbm[1], val=v, type=fbm[2],
                                          reason=attr(meta,"comment"))
            }
        }

        if (x[1] > 0 && recursive) x <- x + autoFilter( verbose=FALSE )
        if (x[1] > 0 && verbose) message(c("Automatic filters have masked",x[1],
              "cells,",x[2],"rows, and",x[3],"cols"), prefix="[-]",
              bgcolor='cyan', color='yellow')
        invisible(x)
    },

    filterByScore = function( min=NA, max=NA, filterEmpty=FALSE, reason=NA ) {
        "\\preformatted{
Apply filters to the current matrix to zero-out cells failing thresholds.
Returns the count of (cells,rows,cols) zeroed (filtered) out.
        min - Minimum allowed value. Cells below this will be set to zero
        max - Maximum allowed value. Cells above it will be set to zero
 filterEmpty - Default FALSE; If true, then the matrix will be 'shrunk' to
              remove rows and columns that are only zeros
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        obj <- matObj()
        rv  <- c(0L, 0L, 0L)
        if (filterEmpty) removeEmpty("Empty rows and cols before score filter")
        if (CatMisc::is.something(min)) {
            ## Zero out entries that fall below min and are not already zero
            fail <- if (min > 0) {
                obj@x < min & obj@x != 0
            } else {
                obj@x < min
            }
            numZ <- sum(fail)
            if (numZ > 0) {
                ## At least some cells were zeroed out
                obj@x[ fail ] <- 0
                testTxt <- paste("score <", min)
                numTxt  <- paste(numZ, "Cells")
                type    <- "Val"
                ijz <- .detailZeroedRowCol( obj, fail, testTxt, reason )
                rv  <- rv + c(numZ, ijz)
                if (filterEmpty) {
                    ## Strip empty rows and columns
                    matrixUse <<- obj
                    if (is.something(reason[1])) {
                        testTxt <- paste(testTxt, reason[1])
                    }
                    removeEmpty(testTxt)
                    obj <- matObj()
                }
            }
            .addAppliedFilter("SCORE", min, reason, '<')
        }
        if (CatMisc::is.something(max)) {
            ## Zero out entries that are above max and not already zero
            fail <- if (max < 0) {
                obj@x > max & obj@x != 0
            } else {
                obj@x > max
            }

            numZ <- sum(fail)
            if (numZ > 0) {
                ## At least some cells were zeroed out
                obj@x[ fail ] <- 0
                testTxt <- paste("score >", max)
                numTxt  <- paste(numZ, "Cells")
                ijz <- .detailZeroedRowCol( obj, fail, testTxt, reason )
                rv  <- rv + c(numZ, ijz)
                if (filterEmpty) {
                    ## Strip empty rows and columns
                    matrixUse <<- obj
                    if (is.something(reason[1])) {
                        testTxt <- paste(testTxt, reason[1])
                    }
                    removeEmpty(testTxt)
                    obj <- matObj()
                }
            }
            .addAppliedFilter("SCORE", max, reason, '>')
        }
        if (rv[1] != 0 && !filterEmpty) matrixUse <<- obj
        invisible(rv)
    },

    .detailZeroedRowCol = function( obj, fail, metric, reason ) {
        "\\preformatted{
Internal method. Takes a SparseMatrix and a (previously calculated)
logical vector of 'failed' rows, and determines which (if any) rows
and columns went from populated to unpopulated (all zeroes). Adds
entries to $filterLog to record them
}"
        ## What rows / columns were affected?
        ij   <- which( fail )
        ## Check if any rows are now totally zeroed out
        i    <- unique( obj@i[ ij ] )
        iBye <- vapply(i, function(xi) all(obj@x[ obj@i == xi ] == 0), TRUE)
        ibs <- sum(iBye)
        ## +1 -> Because Matrix Row/Col are zero-indexed!
        if (ibs > 0) {
            ## Some rows are now all zero because of the filter
            ids <- rownames(obj)[ i[iBye] + 1 ]
            .filterDetails(id=ids, type="Row", metric=metric, reason=reason)
        }
        ## Check the columns
        j    <- unique( obj@j[ ij ] )
        jBye <- vapply(j, function(xj) all(obj@x[ obj@j == xj ] == 0), TRUE)
        jbs  <- sum(jBye)
        if (jbs > 0) {
            ## Some rows are now all zero because of the filter
            ids <- colnames(obj)[ j[jBye] + 1 ]
            .filterDetails(id=ids, type="Col", metric=metric, reason=reason)
        }
        c(ibs, jbs)
    },

    filterByFactorLevel = function( x, keep=TRUE, ignore.case=TRUE,
                                   filterEmpty=FALSE, reason=NA) {
        "\\preformatted{
For factor matrices, filter out values by factor level.
Returns the number of cells zeroed (filtered) out.
          x - A vector of factor levels to remove, either integers or
              character values
       keep - Default TRUE, which will result in keeping only the values
              that match x. If FALSE, the values in x are excluded.
 ignore.case - Default TRUE, which allows factor levels to be match
              while disregarding capitalization
 filterEmpty - Default FALSE; If true, then the matrix will be 'shrunk' to
              remove rows and columns that are only zeros
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        if (!is.factor()) {
            err("Can not filterByFactorLevel() - matrix is not a factor")
            return(invisible(NA))
        }
        obj  <- matObj()
        rv   <- c(0L, 0L, 0L)
        xVal <- integer()
        lvls <- levels()
        nlvl <- length(lvls)
        if (is.numeric(x)) {
            ## Levels appear to be specified by index
            xVal <- if (!is.integer( x )) {
                xVal <- as.integer( x )
                if (!all.equal(x, xVal)) {
                    ## Uh. Probably a problem.
                    err(c("Integer conversion in filterByFactorLevel() appears to have altered your request"), fatal=TRUE)
                }
                xVal
            } else {
                x
            }
            xVal <- unique(xVal)
            low  <- xVal < 1L
            if (sum(low)) {
                bad <- unique(xVal[low])
                xVal <- xVal[ !low ]
                err(c("Rejecting low indices from filterByFactorLevel(): ",
                      paste(bad, collapse=',')))
            }
            high <- xVal > nlvl
            if (sum(high)) {
                bad <- unique(xVal[high])
                xVal <- xVal[ !hig ]
                err(c("Rejecting high indices from filterByFactorLevel(): ",
                      paste(bad, collapse=',')))
            }
        } else if (is.character(x)) {
            ## If we ignore case, a single x value might match
            ## multiple level values. Use a logical matrix to handle
            ## both matches and error reporting
            matchMat <- if (ignore.case) {
                ul <- toupper(lvls)
                apply(matrix(toupper(x), nrow=1), 2, function(z) z == ul)
            } else {
                apply(matrix(x, nrow=1), 2, function(z) z == lvls)
            }
            bad <- base::colSums(matchMat) == 0
            if (sum(bad)) {
                ## At least some values of x could not be matched
                err(c("Some input values could not be matched to levels in filterByFactorLevel(): ",
                      paste(bad, collapse=',')))

            }
            ## The specified-and-valid level indicies are then the
            ## rows with at least one match in the matrix:
            xVal <- which(base::rowSums(matchMat) > 0)
        } else {
            err(c(sprintf("filterByFactorLevel() input error: provided %s.",
                          storage.mode(x)), "Expected  character or numeric"))
            return( invisible(NA) )
        }
        if (length(xVal) == 0) {
            err("No valid levels provided to filterByFactorLevel()")
            return( invisible(NA) )
        }

        ## Even if the user asked to keep a set of levels, the set
        ## being excluded might be smaller (and vice versa). Determine
        ## what the most compact representation of the filter is
        k      <- keep
        lnames <- if (length(xVal) * 2 > nlvl) {
            ## Briefer to refer to the un-specified set
            k <- !k
            lvls[ -xVal ]
        } else {
            lvls[ xVal ]
        }
        op <- if (k) {"=="} else {"!="}
        ## The filter text describes failing cells
        testTxt <- sprintf("Level %s %s", op,
                           paste(lnames, collapse=', '))
        .addAppliedFilter("LEVELS", lnames, reason, op)
        
        ## Working with the @x value vector for the sparse
        ## matrix. Determine which values are 'unwanted'
        fail <- is.element(obj@x, xVal) & obj@x != 0
        if (keep) fail <- !fail
        numZ <- sum( fail )
        if (numZ > 0) {
            # if (filterEmpty) removeEmpty("Empty rows and cols prior to factor filter")
            ## Zero-out the failing cells
            obj@x[ fail ] <- 0
            numTxt  <- paste(numZ, "Cells")
            type    <- "Val"
            ijz <- .detailZeroedRowCol( obj, fail, testTxt, reason )
            rv  <- rv + c(numZ, ijz)
            if (filterEmpty) {
                ## Strip empty rows and columns
                matrixUse <<- obj
                if (is.something(reason[1])) {
                    testTxt <- paste(testTxt, reason[1])
                }
                removeEmpty(testTxt)
                obj <- matObj()
            }
        }
        if (rv[1] != 0 && !filterEmpty) matrixUse <<- obj
        invisible(rv)
    },

    filterByMetadata = function(key, val, MARGIN=NULL, type="like", op='or',
                                reason=NA) {
        "\\preformatted{
Zero out rows or columns that have metadata matching certain values
Returns the number of cells zeroed (filtered) out.
        key - The metadata key/tag/column to match
        val - The value(s) to match
     MARGIN - Default NULL, which will test matched IDs against both rows
              and columns. Can also be 1, row, 2 or col
       type - Default 'like'. Recognized values are:
                like = Finds matches of word/phrase anywhere in text
                regexp = Per 'like', but will interpret as regular expression
                equal  = Finds only full exact matches
              In addition, 'case' can be added to any type to force
              case-sensitive matches.
         op - Default 'or', which will count a match if any of the values
              Alternative is 'and', which requires all values to match.
              'any' and 'all' are also acceptible.
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        ## Normalize the requested key
        key <- key[1]
        x   <- is.element(tolower(metadata_keys()), tolower(key))
        x   <- which(x)
        if (length(x) == 0) {
            err(paste("$filterByMetadata(key='", key,"') does not match any keys"))
            return(NULL)
        }
        if (length(x) > 1) {
            ## Try case insensitive?
            x   <- is.element(metadata_keys(), key)
            x   <- which(x)
            if (length(x) != 1) {
                err(c("filterByMetadata(key =", key,") matches multiple columns"))
                return(NULL)
            }
        }
        colN <- metadata_keys()[x] # The column name (key) we will use
        metric <- colN
        
        if (!is.null(MARGIN)) {
            ## Normalize the margin
            MARGIN <- MARGIN[1]
            if (grepl('(1|row)', MARGIN, ignore.case=TRUE)) {
                MARGIN <- 'row'
            } else if (grepl('(2|col)', MARGIN, ignore.case=TRUE)) {
                MARGIN <- 'col'
            } else {
                err("$filterByMetadata(MARGIN) should be one of NULL,1,2,row,col",
                    fatal=TRUE)
            }
        }

        type <- type[1]
        ## Ignore case flag - if type includes 'case', then do case sensitive
        ic   <- ifelse(grepl('case', type, ignore.case=TRUE), FALSE, TRUE)

        ## Normalize the matching type mat will be a matrix where the
        ## rows correspond to the rows of the matrixMD metadata table,
        ## and the columns are the value(s) that were passed to this
        ## method.
        mat <- if (grepl('(like|sim)', type, ignore.case=TRUE)) {
            ## Flanking wildcard match, no regular expression
            type <- "LIKE"
            if (ic) {
                sapply(val, function(v) grepl(tolower(v),
                    tolower(matrixMD[[colN]]), fixed=TRUE), simplify='matrix')
            } else {
                sapply(val, function(v) grepl(v, matrixMD[[colN]], fixed=TRUE),
                       simplify='matrix')
            }
        } else if (grepl('(regexp)', type, ignore.case=TRUE)) {
            ## Regular expression
            type <- "RegExp"
            sapply(val, function(v) grepl(v, matrixMD[[colN]],
                                          ignore.case=ic), simplify='matrix')
         } else if (grepl('(equal|exact)', type, ignore.case=TRUE)) {
            type <- "Equals"
            if (ic) {
                sapply(val, function(v) tolower(v) == tolower(matrixMD[[colN]]),
                       simplify='matrix')
            } else {
                sapply(val, function(v) v == matrixMD[[colN]],
                       simplify='matrix')
            }
        } else {
            err(c("$filterByMetadata(type =",type,") is not recognized"),
                fatal=TRUE)
        }
        metric <- c(metric, type)
        
        ## Normalize the operator
        fail <- if (grepl('(or|any)', op, ignore.case=TRUE)) {
            ## If any of the values matched, we count it as a fail
            metric <- c(metric, paste(val, collapse=' | '))
            op <- 'or'
            apply(mat, 1, any)
        } else if (grepl('(and|all)', op, ignore.case=TRUE)) {
            ## All must match to count
            metric <- c(metric, paste(val, collapse=' & '))
            op <- 'and'
            apply(mat, 1, any)
        } else {
            err("$filterByMetadata(op) should be one of or,any,and,all",
                fatal=TRUE)
        }

        rv     <- c(0L, 0L, 0L)
        metric <- paste(metric, collapse=' ')
        ## No matches
        if (sum(fail) == 0) return( rv )

        ## Need to consider the hits one by one to check if they've
        ## already been zeroed out or not, and to collect stats. This
        ## could admittedly be slow for large matrices with many hits
        ## and/or rows/cols.
        obj  <- matObj()
        hits <- matrixMD[["id"]][ fail ]
        if (is.null(MARGIN) || MARGIN == 'row') {
            ## Check for matches in rownames
            rn <- match(hits, rownames(obj))
            rn <- rn[ !is.na(rn) ]
            if (length(rn) > 0) {
                inds <- rn - 1 # Zero-index SparseMatrix
                ## Which of these names has non-zero cells in the matrix?
                ## Some may not due to prior filters
                zeroed <- sapply(inds, function(i) obj@i == i & obj@x != 0,
                                 simplify='matrix')
                zHits  <- base::colSums( zeroed )
                inds   <- inds[ zHits != 0 ]
                if (length(inds) != 0) {
                    ## This operation has impacted one or more rows
                    rv[2] <- length(inds)
                    cells <- base::rowSums(zeroed) != 0
                    rv[1] <- rv[1] + sum(cells)
                    obj@x[ cells ] <- 0
                    .filterDetails(id=colnames(obj)[inds + 1], type="Row",
                                  metric=metric, reason=reason)
                }
            }
        }
        if (is.null(MARGIN) || MARGIN == 'col') {
            ## Check for matches in colnames
            cn <- match(hits, colnames(obj))
            cn <- cn[ !is.na(cn) ]
            if (length(cn) > 0) {
                inds <- cn - 1 # Zero-index SparseMatrix
                ## Which of these names has non-zero cells in the matrix?
                ## Some may not due to prior filters
                zeroed <- sapply(inds, function(i) obj@j == i & obj@x != 0,
                                 simplify='matrix')
                zHits  <- base::colSums( zeroed )
                inds   <- inds[ zHits != 0 ]
                if (length(inds) != 0) {
                    ## This operation has impacted one or more columns
                    rv[3] <- length(inds)
                    cells <- base::rowSums(zeroed) != 0
                    rv[1] <- rv[1] + sum(cells)
                    obj@x[ cells ] <- 0
                    .filterDetails(id=colnames(obj)[inds + 1], type="Col",
                                  metric=metric, reason=reason)
                }
            }
        }
        margTxt <- if(is.null(MARGIN)) { "any" } else { MARGIN }
        .addAppliedFilter("METADATA", val, reason,sprintf("'%s' %s", colN,type),
                          sprintf("[%s %s%s]", margTxt, op,
                                  ifelse(ic, ' ignore.case', '')))


        if (rv[1] != 0) matrixUse <<- obj # Update matrix if changed
        invisible(rv)
    },

    nnZero = function(obj=NULL, ...) {
        "\\preformatted{
Return the count of non-zero elements in the matrix. This is essentially
the number of 'live connections' currently stored.
cell
        obj - Default NULL, which will recover the matrix from matObj(),
              passing ... as well (so you can select the raw matrix if
              desired)
}"
        if (is.null(obj)) obj <- matObj(...)
        Matrix::nnzero( obj )
    },

    rCounts = function(obj=NULL, ...) {
        "\\preformatted{
Return an integer vector of non-zero columns for each row
        obj - Default NULL, which will recover the matrix from matObj(),
              passing ... as well (so you can select the raw matrix if
              desired)
}"
        if (is.null(obj)) obj <- matObj(...)
        Matrix::rowSums(obj != 0)
    },

    populatedRows = function(obj=NULL, ...) {
        "\\preformatted{
Return a logical vector indicating which rows have at least one non-zero
cell
        obj - Default NULL, which will recover the matrix from matObj(),
              passing ... as well (so you can select the raw matrix if
              desired)
}"
        rCounts( obj=obj, ...) != 0
    },

    removeEmptyRows = function(reason=NA) {
        "\\preformatted{
Remove all empty rows (those that only contain zeros). Invisibly returns a
vector of removed IDs.
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        obj     <- matObj()
        isEmpty <- !populatedRows()
        ## populatedRows() returns a named logical vector. Take names
        ## from there
        toss    <- names(isEmpty)[ isEmpty ]
        if (length(toss) != 0) matrixUse <<- obj[ !isEmpty, , drop=FALSE]
        invisible(toss)
    },

    cCounts = function(obj=NULL, ...) {
        "\\preformatted{
Return an integer vector of non-zero rows for each column
        obj - Default NULL, which will recover the matrix from matObj(),
              passing ... as well (so you can select the raw matrix if
              desired)
}"
        if (is.null(obj)) obj <- matObj(...)
        Matrix::colSums(obj != 0)
    },

    populatedCols = function(obj=NULL, ...) {
        "\\preformatted{
Return a logical vector indicating which columns have at least one non-zero
 cell
        obj - Default NULL, which will recover the matrix from matObj(),
              passing ... as well (so you can select the raw matrix if
              desired)
}"
        cCounts( obj=obj, ...) != 0
    },

    removeEmptyCols = function(reason=NA) {
        "\\preformatted{
Remove all empty columns (those that only contain zeros). Invisibly
returns a vector of removed IDs.
     reason - Default NA; If specified, a text value that will be added
              to the $filterLog under the 'reason' column
}"
        obj     <- matObj()
        isEmpty <- !populatedCols()
        ## populatedCols() returns a named logical vector. Take
        ## names from there
        toss    <- names(isEmpty)[ isEmpty ]
        if (length(toss) != 0) matrixUse <<- obj[ , !isEmpty, drop=FALSE]
        invisible(toss)
    },

    removeEmpty = function(reason=NA) {
        "\\preformatted{
Remove all empty rows and columns. Invisibly returns a vector of removed IDs.
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        toss <- removeEmptyRows(reason)
        toss <- c(toss, removeEmptyCols(reason))
        invisible(toss)
    },

    map = function(input=NULL, via=NULL, format="data.frame",
                   ignore.case=TRUE, keep.best=FALSE,
                   column.func=max,
                   collapse=NULL, collapse.name=NULL, collapse.token=',',
                   collapse.score=NULL,
                   collapse.factor=NULL, integer.factor=FALSE,
                   add.metadata=TRUE, warn=TRUE,
                   append.to=NULL, append.col=1L, help=FALSE
                   ) {
        "\\preformatted{
Provide a list of IDs, and map/pivot it from one dimension of the matrix
to the other, following 'connections' defined by non-zero cells. Returns
a data.frame with Input and Output columns, plus Score and/or Factor columns.

If collapse='out', rownames will be set to the Output column. Otherwise,
if an input ID results in a unique row, the rowname will be that input ID.
Other-otherwise, the rowname will be the input ID plus '#1', '#2' etc for
each row with that ID. Non-unique input IDs will never generate an
unaltered rowname.

      input - A vector of 'input' IDs. Required, unless append.to is set, in
              which case it will be taken from that table using append.col
        via - Specify if the input matches the 'rows' or 'columns' of the
              matrix. If NULL (default) then your input will be compared to
              the row and column names, and the one with the most matches
              will be chosen (defaulting to 'row' in the event of equal matches)
     format - Default 'data.frame', specifies the output format. Can also be:
                 Vector: Will return a named character vector, with values
                         corresponding to $Output and names as $Input
 ignore.case - Default TRUE, which ignores the capitilazation of IDs
  keep.best - Default FALSE. If TRUE, then only the top-scored cell(s) will
              be kept
 column.func - Default max. If ignore.case is true, it is possible that an
              input ID can match multiple matrix IDs. In this case, multiple
              matching rows will be returned for one ID. column.func is
              applied to reduce this to a single row.
   collapse - Default FALSE, which will cause every pairwise connection to be
              reported. If 'in', then the Input column of the data.frame will be
              unique - any input value that results in multiple output values
              will result in the Output IDs and Score being 'collapsed' to a
              single value (see the collapse.* options below). A value of 'out'
              will do the same, but causes the Output column to be unique, and
              Input and Score are collapsed.
 collapse.name - Default NULL, which will cause multiple names to be
              concatenated into a single value using paste(). Alternatively,
              a user function can be provided. This package also includes the
              crude utility function takeLowestThing() (see documentation).
 collapse.token - Default ',', a string used to concatenate collapsed IDs
              when using paste (collapse.name=NULL)
 collapse.score - Default NULL. Optional function to apply to the $Score
              column when two or more rows are being collapsed into one. The
              function should take a numeric vector as input and return a
              single numeric value. If NULL, the function will be mean(),
              unless the matrix is a pseudo-factor, in which case the object
              method $.autoLevel() will be used. This will generate new
              'hybrid' factors as needed. collapse.token will be used as the
              string to join factors into a new level name.
 integer.factor - Default FALSE, which will cause the Score column to be absent
              and a Factor column (with level values as characters) to be
              present instead. If TRUE then ONLY a Score column (representing
              integer values, perhaps including new likely-meaningless hybrid
              values from $.autoLevel()) will be added. If NULL, then both
              Score and Factor columns are present.
 add.metadata - Default TRUE, which will add all metadata columns that have
              at least some information. FALSE will prevent adding metadata,
              and a character vector will add those specific columns (which
              is up to the user to confirm they exist in the metadata store)
       warn - Default TRUE, which will show warning text if matches failed to
              be made for the input. This information is also always captured
              in attributes attached to the returned data.frame
  append.to - Default NULL. If a data.frame-compliant object, it will become
              the return value, with the mapped columns being added on.
              If collapse='out' then the Output column will be used to
              merge, otherwise the Input column will be utilized. The merge
              column from the provided data.frame is set with append.col
 append.col - The column to use in append.to for merging. Default is 1,
              can provide another column number or name.
}"

        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        ## If we are appending to another data.frame, always take the
        ## input as the designated pivot column:
        if (!is.null(append.to)) input <- append.to[[ append.col ]]
        if (!is.vector(input) || !is.character(input[1])) {
            err("$map() must be provided with a character vector as input")
            return(NA)
        }
        obj   <- matObj()
        ## Build some named vectors to generalize mapping IDs between
        ## input and the matrix
        rn    <- rownames(obj)
        cn    <- colnames(obj)
        inp   <- input
        if (ignore.case) {
            rn  <- setNames(tolower(rn), rn)
            cn  <- setNames(tolower(cn), cn)
            inp <- setNames(tolower(input), input)
        } else {
            names(rn)  <- rn
            names(cn)  <- cn
            names(inp) <- input
        }

        ## Is the input unique, or duplicated?
        dupIn <- base::duplicated(inp)
        if (sum(dupIn) == 0) {
            dupIn <- character()
        } else {
            ## Some of the input is duplicated
            dupIn <- unique( inp[dupIn] )
            dupIn <- vapply(dupIn, function(x) {
                paste(names(rn)[which(rn == x)], collapse=collapse.token)
            }, "")
            inp  <- inp[ !base::duplicated(inp) ]
        }

        ## Work out if we're collapsing anything
        colIn  <- FALSE
        colOut <- FALSE
        if (!is.null(collapse)) {
            ## Assure the collapse token is a single string:
            collapse.token <- collapse.token[1]
            if (grepl('in', collapse[1], ignore.case=TRUE)) {
                colIn  <- TRUE
            } else if (grepl('out', collapse[1], ignore.case=TRUE)) {
                colOut <- TRUE
            } else {
                err("mapToCol() collapse parameter must be 'in' or 'out'")
            }
        }

        ## Determine which dimension of the matrix we're matching to
        if (is.null(via)) {
            ## Automatically determine if we are starting with rows/cols
            rCnt <- sum(is.element(inp, rn))
            cCnt <- sum(is.element(inp, cn))
            if (cCnt > rCnt) {
                ## We have more matches in columns than rows
                via   <- "col"
            } else {
                via   <- 'row'
            }
        } else if (grepl('^row', via[1], ignore.case=TRUE)) {
            via   <- 'row'
        } else if (grepl('^col', via[1], ignore.case=TRUE)) {
            via   <- 'col'
        } else {
            err(c("mapToCol() unrecognized 'via' argument.",
                  "Should be 'row' or 'col', or NULL for automatic"))
            return(NA)
        }

        ## What are the relavent IDs we will be comparing to?
        matIds <- if (via == 'col') {
            ## We're entering the matrix from columns
            obj <- matObj(transpose=TRUE) # transpose for simplicity below
            cn
        } else {
            rn # Matching against rows
        }
        ## What are the original names (prior to any to.lower() operations):
        matNms <- names(matIds)
        inpNms <- names(inp)
        
        ## Determine any unknown IDs
        unk <- base::setdiff(inp, matIds)
        ids <- if (length(unk) > 0) {
            ## Some user IDs do not have a match in the matrix.
            unk <- inpNms[ match(unk, inp) ] # Restore user's case:
            base::intersect(inp, matIds)
        } else {
            inp
        }
        
        unMap   <- character() # Unmapped IDs (known, but no path to target)
        dupMat  <- character() # Matrix IDs in two or more rows (after to.lower)
        multIn  <- character() # Input terms that have multiple Outputs
        isFac   <- is.factor() # Just a handy boolean

        ## Build the data.frame by columns
        vecStp <- 1000              # Vector chunking amount
        inpCol <- character(vecStp) #   Vector of input IDs
        outCol <- character(vecStp) #   Vector of output IDs
        scrCol <- numeric(vecStp)   #   Vector of scores
        rcnt   <- 0                 #   Number of rows in vectors so far

        ## Wrap up the function used to collapse scores here, in order
        ## to collect all the tests in one place and pre-choose which
        ## function to use. Each function will take a single vector of
        ## numeric values
        valueCollapseFunction <- if (!is.null(collapse.score)) {
            ## The user is supplying their own function to handle
            ## value (score) collapse
            collapse.score
        } else if (isFac) {
            ## Default method to handle factors by making
            ## new 'hybrid' ones
            function (vec) .autoLevel( vec, sep=collapse.token )
        } else {
            mean
        }
        ## Likewise, the function used to collapse names/IDs when two
        ## or more rows need to be stuffed into one. Also takes a
        ## single vector, in this case of characters.
        nameCollapseFunction <- if (is.null(collapse.name)) {
           # If no user-supplied function is given, then just use
           # paste() with collapse.token as the separator.
           function(vec) paste(vec, collapse=collapse.token)
        } else {
            ## User-supplied function
            collapse.name
        }


        ## Cycle through each input id, building columns as we go:
        for (id in ids) {
            idNm <- inpNms[ match(id, inp) ] # User's name (before to.lower)
            # Select the indices for matrix row(s) matching the id:
            ind  <- which(id == matIds)
            rows <- obj[ ind, , drop=FALSE]
            ## Select the columns that are non-zero
            cSum <- Matrix::colSums( rows )
            rows <- rows[ , cSum != 0, drop=FALSE]
            if (length(rows) == 0) {
                ## No mappings! The input ID existed in the matrix,
                ## but there are no non-zero mappings to the other
                ## edge.
                rows  <- matrix(0L, 1, 1, dimnames=list(r=NA, c=NA))
                unMap <- c(unMap, idNm)
                ## next
            }
            ## Move from a subset of the matrix to a single vector
            vec <- if (nrow(rows) > 1) {
                ## We are matching multiple matrix rows with this
                ## ID. This can happen when there are multiple cases
                ## for the same rowname (eg gene symbols "p40" and
                ## "P40")
                dupMat <- c(dupMat,paste(rownames(rows),
                                         collapse=collapse.token))
                apply(rows, 2, column.func)
            } else {
                ## The '[' opperator is not honoring column names in
                ## Matrix objects so I need to explicitly setNames:
                setNames(rows[1, , drop=TRUE], colnames(rows))
                ## This is maybe because dimensions in dgTMatrix
                ## objects are stored in slots, not attributes?
            }

            if (keep.best) {
                ## Only keep the top-scored result(s)
                mx  <- max(vec)
                vec <- vec[ vec == mx ]
            }

            ## Do we end up with more than one output term?
            if (length(vec) != 1) {
                multIn <- c(multIn, idNm)
                if (colIn) {
                    ## Request to collapse by input ID. We can do that
                    ## here, since ids represents all unique input IDs
                    nmVal  <- nameCollapseFunction( names(vec) )
                    colVal <- valueCollapseFunction(vec)
                    vec    <- setNames(colVal, nmVal)
                }
            }
            vl <- length(vec) # How many results we have
            need <- rcnt + vl - length(inpCol) # Do we have enough space?
            if (need > 0) {
                ## No, we need to grow our results vectors
                inpCol  <- c(inpCol, character(need + vecStp))
                outCol  <- c(outCol, character(need + vecStp))
                scrCol  <- c(scrCol, numeric(need + vecStp))
            }
            newInds <- seq(rcnt+1, length.out=vl) # Indices to 'inject' into
            inpCol[ newInds ] <- rep(idNm, vl)    #   Input Id, replicated
            outCol[ newInds ] <- names(vec)       #   Output IDs
            scrCol[ newInds ] <- vec              #   Scores
            rcnt <- rcnt + vl                     #   Current number of rows
        }
        if (length(unk) > 0) {
            ## Push the unknown IDs onto the end
            vl <- length(unk)
            need <- rcnt + vl - length(inpCol) # Do we have enough space?
            if (need > 0) {
                ## No, we need to grow our results vectors
                inpCol  <- c(inpCol, character(need))
                outCol  <- c(outCol, character(need))
                scrCol  <- c(scrCol, numeric(need))
            }
            newInds <- seq(rcnt+1, length.out=vl) # Indices to 'inject' into
            inpCol[ newInds ] <- unk         # Each unknown ID
            outCol[ newInds ] <- rep(NA, vl) #   Output IDs = NA
            scrCol[ newInds ] <- rep(NA, vl) #   Scores = NA
            rcnt <- rcnt + vl                #   Current number of rows
        }
        sl <- seq_len(rcnt)
        rv <- data.frame(Input  = inpCol[sl],
                         Output = outCol[sl],
                         Score  = scrCol[sl],
                         stringsAsFactors=FALSE)
        if (!isFac || all(integer.factor)) {
            ## Add a numeric score column
            if (isFac) {
                rv$Score <- as.integer(scrCol[sl])
            } else {
                rv$Score <- scrCol[sl]
            }
        }

        multOut <- base::duplicated(rv$Output)
        if (sum(multOut) == 0) {
            multOut <- character()
        } else {
            ## Some of the output terms come from multiple inputs
            inds    <- which(multOut)
            ## Note them so they can be reported to user:
            multOut <- unique( rv$Output[multOut] )
            if (colOut) {
                ## Request to collapse by input ID
                for (nm in multOut) {
                    ## Get all rows for this Output name:
                    nInds <- which(rv$Output == nm)
                    ## We will keep the first row:
                    keep  <- nInds[1]
                    ## Aggregate the scores, stuff into the first index:
                    rv[keep, "Score"] <- valueCollapseFunction(scrCol[ nInds ])
                    ## ... and the input names:
                    rv[keep,"Input"]  <- nameCollapseFunction( inpCol[ nInds ])
                }
                ## Delete the "extra" rows
                rv   <- rv[ -inds, ]
                ## Update length information:
                rcnt <- nrow(rv)
                sl   <- seq_len(rcnt)

            }
        }

        ## Set rownames
        rownames(rv) <- if (colOut) {
            ## Collapsing on output -> Name rows by Output column
            rv$Output
        } else if (colIn || length(multIn) == 0) {
            ## Collapsing on Input, or Input (by chance or design) is unique:
            rv$Input
        } else {
            strict.unique( rv$Input )
        }
        
        if (isFac && !any(integer.factor)) {
            ## Add factor text column
            inds <- as.integer(rv$Score)
            ## If we recognized the input but failed to find a map,
            ## the score is zero. This will evaluate as character(0)
            ## inside lvlVal[], so map these values to NA instead:
            inds[ inds == 0 ] <- NA
            if (all(is.na(inds))) {
                ## A single NA can be a bit of a pain, manually set:
                rv$Factor <- as.character( rep(NA, length(inds)) )
            } else {
                rv$Factor <- lvlVal[ inds ]
            }
        }

        if (CatMisc::is.something(add.metadata)) {


### TODO: Need to manage cases where colIn is true, and Output is a
### token-separated string of IDs. That is, need to grab metadata in
### some fashion for these synthetic values. Should probably hold a
### temp list for each row of original values
            
            ## Request to include metadata columns
            ## Map our output IDs to the Metadata data.table :
            md  <- matrixMD[rv$Output, setdiff(colnames(matrixMD), "id"),
                            with=FALSE]
            mdc <- colnames(md) # All available metadata columns
            ## If the param is not a character vector, take all
            if (!is.character(add.metadata)) add.metadata <- mdc
            for (col in add.metadata) {
                if (is.element(col, mdc)) {
                    ## This is a known metadata column. Include it,
                    ## unless it's all NAs
                    if (!all(is.na(md[[ col ]]))) rv[[col]] <- md[[ col ]]
                } else {
                    ## Column does not exist, put a blank text column
                    ## in the output
                    rv[[col]] <- character(rcnt)
                }
            }
        }

        ## Also attach a brief explanation of each attribute:
        Notes <- list(
            Via = "Whether your input was matched to rows or columns of the matrix",
            'Dup.In' = "IDs that were present twice or more in input (possibly after case removal)",
            'Dup.Mat' = "IDs that were present twice or more in matrix (after case removal)",
            'Mult.In' = "Input IDs that generated multiple output values",
            'Mult.Out' = "Output IDs that were generated from multiple inputs",
            Unmapped = "Input ID is also in the matrix, but does not have a target with non-zero score. Score will be zero",
            Unknown = "Input ID could not be matched to any in the matrix. Score will be NA")

        if (!is.null(append.to)) {
            ## User wants result columns added to their own data.frame
            mainCol  <- if (colOut) { 1L } else { 2L }
            cnames   <- colnames(rv)
            transfer <- c(mainCol, seq(3, length(cnames)))
            added    <- character()
            for (i in transfer) {
                ## Make sure the appended column is unique as we
                ## extend the data.frame:
                cn  <- cnames[i]
                u   <- make.unique(c(colnames(append.to), cn))
                cnu <- u[ length(u) ]
                append.to[[ cnu ]] <- rv[ input, i ]
                added <- c(added, cnu)
            }
            rv <- append.to
            attr(rv, "Appended") <- added
            Notes[["Appended"]]  <- "Column names that were appended to the data.frame you provided with append.to="
        } else if (grepl('vec', format[1], ignore.case=TRUE)) {
            ## Vector format.
            rv <- setNames(rv$Output, rv$Input)
            ## This can result in non-unique names. Do we want to use
            ## strict.unique() to uniquify names?
            class(rv) <- c("mapResult", class(rv))
        }

        ## Attach some summary attributes:
        attr(rv, "Unmapped")   <- unMap
        attr(rv, "Unknown")    <- unk
        attr(rv, "Via")        <- via
        attr(rv, "Dup.In")     <- dupIn
        attr(rv, "Dup.Mat")    <- dupMat
        attr(rv, "Mult.In")    <- multIn
        attr(rv, "Mult.Out")   <- multOut
        if (isFac) attr(rv, "Levels") <- lvlVal
        attr(rv, "Notes")      <- Notes
            
        if (warn) {
            ## Alert user to any non-unique relationships:
            msg <- character()
            cm  <- function(w,x,cl) colorize(paste(w,":",length(x)), cl)
            if (length(unk) != 0)   msg <- c(msg,cm("Unknown",unk,"magenta"))
            if (length(unMap) != 0) msg <- c(msg,cm("Unmapped",unMap,"red"))
            if (length(dupIn) != 0)  msg <- c(msg,cm("Dup.In",dupIn,"blue"))
            if (length(dupMat) != 0) msg <- c(msg,cm("Dup.Mat",dupMat,"cyan"))
            if (length(multIn) != 0) msg <- c(msg,cm("Mult.In",multIn,"yellow"))
            if (length(multOut) != 0) msg <- c(msg,cm("Mult.Out",multOut,"green"))
            if (length(msg) !=0) base::message("Non-unique events: ",
                                               paste(msg, collapse=', '))
        }
        rv
    },

    .autoLevel = function (vals, sep=',', decreasing=TRUE) {
        "\\preformatted{
Designed to be used when map() collapses factorized scores. If two or more
factor values are being collapsed, the function will collapse the levels into
a new token-separated level
       vals - Required, a vector of integer factor values
        sep - Default ',', the token used to separate levels
 decreasing - Default TRUE, which will sort factor values from largest to
              smallest (in order to put the 'best'-scored factor at the
              front of the list)
}"
        ## Setting this as an object method to help access/modify the
        ## lvlVal field
        if (!is.factor()) return(NA)
        if (is.numeric(vals)) {
            vals <- as.integer(vals)
        } else if (!is.integer(vals)) {
            return(NA)
        }
        lvl <- paste(lvlVal[unique(sort(vals, decreasing=decreasing))],
                     collapse=sep)
        ## If this is the first time seeing a hybrid value, add it to the
        ## factor levels:
        if (!is.element(lvl, lvlVal)) lvlVal <<- c(lvlVal, lvl)
        ## Return the index:
        which(lvl == lvlVal)
    },

    .filterDetails = function ( id=NA, type=NA,  metric=NA, reason=NA) {
        if (is.null(metric)) return(NA)
        if (any(is.na(id))) {
            ## Should not happen! This method should only be called on
            ## defined identifiers
            err(sprintf(
                "Filtering noted for metric='%s' on type='%s', but some IDs are NA",
                if (all(is.null(metric))) { "-NULL-" } else { metric },
                if (all(is.null(type))) { "-NULL-" } else { type }),
                prefix = "[CodeError]")
            id <- id[ !is.na(id) ]
        }
        if (length(id) == 0) {
            ## Should not happen! There should always be at least one
            ## identifier
            err(sprintf(
                "Filtering noted for metric='%s' on type='%s', but no IDs provided",
                if (all(is.null(metric))) { "-NULL-" } else { metric },
                if (all(is.null(type))) { "-NULL-" } else { type }),
                prefix = "[CodeError]")
            return(NA)
        }
        row <- data.table::data.table(id=as.character(id), metric=metric,
                                      type=type, reason=reason, key="id")
        ## Do not add any entries that are already recorded as
        ## filtered (only note the first exclusion)
        newR <- base::setdiff(row[["id"]], filterLog[["id"]])
        filterLog <<- extendDataTable(filterLog, row[newR])
        data.table::setkeyv(filterLog, "id")
        filterLog
    },

    .addAppliedFilter = function (key, val, com=NULL, pre=NULL, pro=NULL) {
        "\\preformatted{
Internal method, stores filter text in $setFilters
      key - A key indicating the filter type, eg 'SCORE' or ''
      val - Text defining the filter parameters
      com - Optional comment (reason) for the filter
      pre - Optional text (eg operator) to go before the values
      pro - Optional text (eg modifiers) to go after the values
}"
        txt <- c(toupper(key))
        if (CatMisc::is.something(pre)) txt <- c(txt, pre)
        txt <- c(txt, paste(val, collapse=sfSep))
        if (CatMisc::is.something(pro)) txt <- c(txt, pro)
        if (CatMisc::is.something(com)) txt <- c(txt, "##", com)
        txt <- paste(txt, collapse=' ')
        if (!is.element(txt, setFilters)) setFilters <<- c(setFilters, txt)
        txt
    },

    filterSummary = function (reason=TRUE) {
        "\\preformatted{
Tallies number of filtered objects (generally cells) by filter criteria
and optionally reason. Returns a data frame
   reason - Default TRUE, which includes the reason text in the count
            grouping
}"

        rv <- if (reason) {
            dplyr::count(filterLog, metric, type, reason)
        } else {
            dplyr::count(filterLog, metric, type)
        }
        rv <- if (nrow(rv) == 0) {
            ## Nothing filtered
            data.frame(metric=character(), reason=character(),
                       Col=integer(), Row=integer())
        } else {
            reshape2::dcast(rv, metric + reason ~ type, sum,
                            value.var = "n",na.rm=TRUE)
        }
        ## Class to utilize pretty-print function
        class(rv) <- c("mapFilter", class(rv))
        attr(rv, "useCol")  <- useColor()
        attr(rv, "filters") <- setFilters
        rv
    },

    appliedFilters = function (new=NULL) {
        "\\preformatted{
Filters applied to matrix so far, represented as text that can be read
both by humans and computationally.
      new - Provide new filters to apply, using the same format.
            Designed to allow a snapshot of the filter state to be
            reapplied to the matrix after a $reset()
}"
        if (!is.null(new)) {
            err("
ToDo: STILL WORKING ON ROUND-TRIP PARSING FILTER TEXT
")
        }
        setFilters
    },

    is.factor = function () {
        "\\preformatted{
TRUE if the matrix is a pseudo-factor (levels have been defined),
otherwise FALSE
}"
        CatMisc::is.something(lvlVal)
    },

    levels = function( asFactor=FALSE ) {
        "\\preformatted{
Returns factor levels, if appropriate. If not, returns NULL
 asFactor - Default FALSE, which will return an ordered character vector of the
            level values (names). If true, a factor will be returned with
            appropriate levels assigned
}"
        if (is.factor()) {
            if (asFactor) { factor(lvlVal) } else { lvlVal }
        } else {
            NULL
        }
    },

    as.gmt = function( obj=NULL, transpose=FALSE, file=NULL, ... ) {
        "\\preformatted{
Converts matrix into a block of GMT-formatted text
        obj - Default NULL, which will recover the matrix from matObj(),
              passing ... as well. Alternatively a Matrix can be provided.
  transpose - Default FALSE, which will presume that the rows are sets. If TRUE,
              then columns will be taken as sets
       file - Default NULL, if defined then the output will be written to that
              path
}"
        ## http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
        if (is.null(obj)) obj <- matObj(transpose=transpose, ...)
        ## Only keep rows (sets) with at least one object:
        hasData  <- populatedRows(obj)
        obj      <- obj[hasData, , drop=FALSE]
        setNames <- rownames(obj) # Rows are sets
        memNames <- colnames(obj) # Columns are the potential members
        descr    <- matrixMD[.(setNames), "Description", with=FALSE]
        ns       <- length(setNames)
        out      <- character(ns)
        for (i in seq_len(ns)) {
            row    <- obj[i, , drop=TRUE]
            out[i] <- sprintf("%s\n", paste(c(setNames[i], descr[i],
                                           memNames[row != 0]), collapse="\t"))
        }
        if (is.null(file)) {
            paste(out, sep='')
        } else {
            invisible(cat(out, sep='', file=file))
            
        }
    },

 

    .readMatrix = function ( format = "", ... ) {
        if (!CatMisc::is.def(file)) err("AnnotatedMatrix objects must define 'file' when created", fatal = TRUE)
        objFile  <- paste(file,'rds', sep = '.')
        rdsOk    <- rdsIsCurrent( file )
        if (is.na(rdsOk) || rdsOk) {
            ## The RDS file is present and can be used
            dateMessage(paste("Reading serialized matrix",
                              colorize(objFile,"white")), prefix = "  ")
            if (is.na(rdsOk)) message("Original data file not available!",
                                      prefix="[CAUTION]", color="yellow")
            rv <- readRDS(objFile)
            fromRDS <<- TRUE
        } else {
            if (!file.exists(file)) err(c("Can not make AnnotatedMatrix",
                                          "File does not exist : ", file),
                                        fatal = TRUE)
            fname   <- file
            colName <- colorize(file,"white")
            is.gz   <- CatMisc::parenRegExp('(.+)\\.gz$', fname)
            if (!is.na(is.gz[1])) fname <- is.gz[1]
            
            if (grepl('(mtx|matrixmarket)', format, ignore.case = TRUE) ||
                grepl('\\.mtx', fname, ignore.case = TRUE)) {
                dateMessage(paste("Reading MatrixMarket file", colName))
                rv <- parse_MatrixMarket_file( file )
            } else if (grepl('(txt|text)', format, ignore.case = TRUE) ||
                       any(dir.exists(file)) ||
                       grepl('\\.(txt|text|list)$',fname, ignore.case = TRUE)) {
                dateMessage(paste("Reading simple list from file", colName))
                rv <- parse_Text_file( file=file, ... )
            } else if (grepl('(lol)', format, ignore.case = TRUE) ||
                       grepl('\\.(inp)$', fname, ignore.case = TRUE)) {
                dateMessage(paste("Reading List-of-Lists", colName))
                rv <- parse_ListOfLists_file( file )
            } else if (grepl('(gmt)', format, ignore.case = TRUE) ||
                       grepl('\\.(gmt)$', fname, ignore.case = TRUE)) {
                dateMessage(paste("Reading lists from GMT file", colName))
                rv <- parse_GMT_file( file=file, ... )
            } else {
                err(c("Can not make AnnotatedMatrix - unrecognized file type: ",
                      file), fatal = TRUE)
            }
            ## Stub DT if no metdata was defined:
            if (!CatMisc::is.something(rv$metadata)) rv$metadata <-
                    data.table::data.table( id = character(), key = "id" )
            dateMessage(paste("Serializing matrix to file",
                              colorize(objFile,"white")), prefix = "  ")
            saveRDS(rv, objFile)
        }
        if (is.null(rv)) err(c("Failed to read file:", file), fatal = TRUE)

        matrixRaw  <<- rv$matrix
        matrixMD   <<- rv$metadata
        if (!is.null(rv$rowChanges)) rowChanges <<- rv$rowChanges
        if (!is.null(rv$colChanges)) colChanges <<- rv$colChanges
        if (CatMisc::is.def(rv$levels)) lvlVal  <<- rv$levels
        ## Set default parameters, without clobbering any already set
        if (CatMisc::is.def(rv$params)) setParamList(rv$params, clobber = FALSE)

        ## Numeric conversion to prevent integer overflow on product
        cnum       <- as.numeric(ncol(matrixRaw))
        rnum       <- as.numeric(nrow(matrixRaw))
        nnz        <- nnzero(matrixRaw)
        pnz        <- smallNumberFormatter( nnz / (cnum * rnum) )
        actionMessage(sprintf("%d x %d matrix, %s non-zero", rnum, cnum, pnz),
                      prefix = "  ")
        rv
    },

     metadata_keys = function ( id = NULL, key = NULL) {
        "\\preformatted{
Returns a character vector of metadata key names (eg 'Description')
}"
        ## Exclude the id column
        base::setdiff(names(matrixMD), "id")
     },

     metadata = function ( id=NULL, key=NULL, na.rm=TRUE, drop=TRUE,
                          verbose=TRUE) {
        "\\preformatted{
Select metadata by id, key or both
         id - Optional vector of IDs to query. If not provided, all IDs
              will be returned. Bear in mind that the metadata holds both
              row and column IDs mixed together
        key - Optional vector of key/tag names, if not provided then all
              available ones will be returned.
      na.rm - Default TRUE, which will remove NA values from returned
              results. Does not apply when both id and key are specified.
       drop - Default TRUE. If FALSE, the return value will be a data.table
              If TRUE, and zero or one metadata columns are present, then
              a named vector is returned.
    verbose - Default TRUE, which will warn about certain issues
}"
        rv <- NULL
        ## All metadata if neither key nor id is specified
        if (is.null(id) && is.null(key)) return( matrixMD )
        if (!is.null(key)) {
            ## data.table is unhappy if asked for non-existent columns
            mkeys   <- metadata_keys()
            key     <- base::setdiff(key, "id") # Initially assure id is not present
            unknown <- base::setdiff(key, mkeys)
            if (length(unknown) > 0) {
                key <- intersect(key, mkeys)
                if (verbose) err(c("Request for unknown metadata key(s):",
                                   unknown))
                if (length(key) == 0) {
                    ## No more keys survive. Are we done here?
                    key <- NULL
                    if (is.null(id)) return( rv )
                }
            }
            ## Make sure the id column is included
            useKey <- c("id", key)
        }
        ## data.table does *not* like numeric keys... They're treated
        ## as (fluid) row numbers
        if (!is.null(id)) id <- as.character(id)

        if (is.null(id)) {
            ## Key request only
            rv <- matrixMD[ , useKey, with=FALSE ]
            if (na.rm) {
                ## Remove null rows. Find NA values in columns:
                naCol <- sapply(useKey, function(x) is.na( rv[[ x ]] ) )
                ## And now find rows that are all() NAs:
                naRow <- apply(naCol, 1, all)
                if (any(naRow)) {
                    ## ... and if any columns are all NAs, remove them:
                    rv <- rv[ !naRow, ]
                    ## Ok to use logicals when subsetting data.table *rows*
                }
            }
        } else if (is.null(key)) {
            ## ID request only. Do not use with=FALSE !
            rv <- matrixMD[ id,  ]
            if (na.rm) {
                ## Remove null columns. Find NA values in columns:
                cols   <- colnames(rv)
                naCol1 <- sapply(cols, function(x) is.na( rv[[ x ]] ) )
                ## And then columns that are all null:
                naCol2 <- apply(naCol1, 2, all)
                if (any(naCol2)) {
                    ## ... and if any columns are all NAs, remove them:
                    rv <- rv[ , which(!naCol2), with=FALSE ]
                    ## Must use column indices or names, and must use with=F !
                }
            }
        } else {
            ## Both
            rv <- matrixMD[ id, useKey, with=FALSE ]
            ## We will not remove either rows or columns, since both
            ## were explicitly requested
        }

        cns <- colnames(rv)
        if (drop && length(cns) <= 2) {
            names <- rv[[ "id" ]]
            if (length(cns) == 1) {
                ## Big batch of nothing (no metadata cols)
                rv <- setNames(rep(NA, length(names)), names)
            } else {
                ## What metadata column do we have?
                mCol <- setdiff(cns, "id")
                rv   <- setNames(rv[[ mCol[1] ]], names)
            }
        }
        rv
    },

    getRow = function ( x = NA, format = "vector", sort = NA,
        min = NA, max = NA, ...) {
        obj <- matObj(...)
        fmt <- tolower(substr(format, 1, 1))
        ## Matrix:: does not handle unrecognized names gracefully
        invalid <- NULL
        if (is.character(x)) {
            valid <- is.element(x, rownames(obj))
            if (sum(valid) < length(x)) {
                invalid <- x[ !valid ]
                x <- x[ valid ]
            }
        }
        if (length(x) == 0) {
            ## No valid selectors
            rv <- if (fmt == "v") {
                ## Empty vector
                c()
            } else if (fmt == "d") {
                if (is.null(invalid)) {
                    data.frame() # plain empty DF
                } else {
                    ## Data frame with columns named for invalid requests
                    data.frame(sapply(invalid, function(x) NA, simplify=FALSE))
                }
            } else {
                if (is.null(invalid)) {
                    matrix() # plain empty matrix
                } else {
                    matrix(numeric(), ncol=length(invalid),
                           dimnames=(list(c(), invalid)))
                }
            }
            return(rv)
        }

        ## Get the rows we are interested in:
        rv  <- obj[ x, , drop=FALSE]
        ## Apply (by "masking" to zero) min/max limits, if requested
        if (!is.na(min)) rv[ rv < min ] <- 0
        if (!is.na(max)) rv[ rv > max ] <- 0
        ## Eliminate all zero cols: https://stackoverflow.com/a/6632287
        rv  <- rv[ , Matrix::colSums(abs(rv)) != 0, drop=FALSE ]
        
        
        if (fmt == "v") {
            ## vector
            ## Just return a vector of matching column names
            colnames(rv)
        } else if (fmt == "d") {
            ## data.frame
            df <- as.data.frame(t(as.matrix(rv)))
            df$colName = colnames(rv)
            df
        } else {
            ## Just return the matrix, transposed
            Matrix::t(rv)
        }
    },

    getCol = function ( x = NA, format = "vector", sort = NA,
        min = NA, max = NA, ...) {
        obj <- matObj(...)
        fmt <- tolower(substr(format, 1, 1))
        ## Matrix:: does not handle unrecognized names gracefully
        invalid <- NULL
        if (is.character(x)) {
            valid <- is.element(x, colnames(obj))
            if (sum(valid) < length(x)) {
                invalid <- x[ !valid ]
                x <- x[ valid ]
            }
        }
         if (length(x) == 0) {
            ## No valid selectors
            rv <- if (fmt == "v") {
                ## Empty vector
                c()
            } else if (fmt == "d") {
                if (is.null(invalid)) {
                    data.frame() # plain empty DF
                } else {
                    ## Data frame with columns named for invalid requests
                    data.frame(sapply(invalid, function(x) NA, simplify=FALSE))
                }
            } else {
                if (is.null(invalid)) {
                    matrix() # plain empty matrix
                } else {
                    matrix(numeric(), ncol=length(invalid),
                           dimnames=(list(c(), invalid)))
                }
            }
            return(rv)
        }
        ## Get the columns we are interested in:
        rv  <- obj[ , x, drop=FALSE ]
        ## Apply (by "masking" to zero) min/max limits, if requested
        if (!is.na(min)) rv[ rv < min ] <- 0
        if (!is.na(max)) rv[ rv > max ] <- 0
        ## Eliminate rows of all zeros: https://stackoverflow.com/a/6632287
        rv  <- rv[ Matrix::rowSums(abs(rv)) != 0, , drop=FALSE ]
        
        if (fmt == "v") {
            ## vector
            ## Just return a vector of matching row names
            rownames(rv)
        } else if (fmt == "d") {
            ## data.frame
            df <- as.data.frame(as.matrix(rv))
            df$colName = rownames(rv)
            df
        } else {
            ## Just return the matrix
            rv
        }
    },

    show = function (...) { cat( .self$matrixText(...) ) },

    matrixText = function ( pad = "", useObj=NULL, fallbackVar=NULL,
                          compact=FALSE, color=NULL ) {
        ## Check for stub object created when calling help on the base class:
        if (CatMisc::is.empty.field(EvLogObj)) return("")

        if (is.null(color)) color <- useColor() # Use EventLogger setting
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        ## Variable name for use in sample methods
        objName <- .self$.selfVarName("myMatrix", fallbackVar)
        whtName <- doCol(objName, "white")
        msg <- doCol(sprintf("%s sparse matrix\n", class(.self)[1]), "blue")
        ## Even if we read the file from RDS, we will base the
        ## reported age on the original data file (if available)
        src <- ifelse(file.exists(file), file, sprintf("%s.rds", file))
        age <- sprintf("%.1f",difftime(Sys.time(),file.mtime(src),units="days"))
        msg <- sprintf("%s    File: %s%s [%s days old]\n",
                       msg, doCol(file, "white"),
                       ifelse(fromRDS, doCol('.rds', "purple"),""),
                       doCol(age,"red"))
        name <- param("name")
        if (CatMisc::is.something(name))
            msg <- sprintf("%s    Name: %s\n", msg, doCol(name, "white"))
        if (!compact) {
            src <- param("source")
            if (CatMisc::is.something(src))
                msg <- paste(msg, sprintf("  Source: %s\n", doCol(src,"white")),
                             collapse='', sep='')
            
            desc <- param("description")
            if (CatMisc::is.something(desc))
                msg <- paste(msg, sprintf("    Desc: %s\n",doCol(desc,"white")),
                             collapse='', sep='')
            
            sd <- param("scoredesc")
            if (CatMisc::is.something(sd))
                msg <- paste(msg, sprintf("  Scores: %s\n", doCol(sd, "white")),
                             collapse='', sep='')
            
            if (is.factor()) {
                ## Report the factor levels
                indent <- "\n      ";
                lvl <- paste(doCol(strwrap(paste(lvlVal, collapse = ', '),
                                           width = 0.7 * getOption("width")),
                                   "purple"), collapse = indent)
                msg <- sprintf("%s    %s%s%s\n", msg, doCol(
                      "Values are factors:", "yellow"), indent, lvl)
            }
        }
        
        ## If an external "used" matrix is not provided, use internal field:
        mat  <- if (is.null(useObj)) { matObj() } else { useObj }
        dimNames <- names(dimnames(mat))
        popr <- populatedRows(mat)
        nr   <- as.numeric( sum(popr) )
        popc <- populatedCols(mat)
        nc   <- as.numeric( sum(popc) )
        nz   <- nnZero(mat)

        ## Note number, namespace and examples for rows:
        rdn  <- dimNames[1]
        rN   <- rownames(mat)[popr]
        rblk <- sprintf("  %8d", nr)
        if (CatMisc::is.something(rdn)) rblk <- paste(rblk, doCol(rdn, "cyan"))
        rblk <- paste(rblk, "rows")
        unpr <- sum(!popr)
        if (unpr > 0) rblk <- paste(rblk, sprintf("(+%d empty)", unpr))
        rblk <- paste(rblk, "eg:",
                      doCol(substr(paste(rN[1:pmin(3,length(rN))],
                                         collapse=', '), 1, 60), "cyan"))
        msg <- paste(msg, rblk, "\n", sep='')

        ## Note same for columns:
        cdn  <- dimNames[2]
        cN   <- colnames(mat)[popc]
        cblk <- sprintf("  %8d", nc)
        if (CatMisc::is.something(cdn)) cblk <- paste(cblk, doCol(cdn, "cyan"))
        cblk <- paste(cblk, "cols")
        unpc <- sum(!popc)
        if (unpc > 0) cblk <- paste(cblk, sprintf("(+%d empty)", unpc))
        cblk <- paste(cblk, "eg:",
                      doCol(substr(paste(cN[1:pmin(3,length(cN))],
                                         collapse=', '), 1, 60), "cyan"))
        msg <- paste(msg, cblk, "\n", sep='')

        ## Note number of non-zero cells
        perPop <- smallNumberFormatter( nz / (nr * nc))
        msg <- paste(msg, sprintf("  %d non-zero cells (%s)\n", nz,
                                  doCol(perPop, "red")))

        ## Note if we had to remap one or more R/C names
        changed <- character()
        if (CatMisc::is.def(rowChanges)) changed <- c(changed, doCol(
            paste(length(rowChanges), "Rows ($rowChanges)"), "yellow"))
        if (CatMisc::is.def(colChanges)) changed <- c(changed, doCol(
            paste(length(colChanges), "Cols ($colChanges)"), "yellow"))

        if (length(changed) != 0) msg <- sprintf(
                      "%s    %s names needed to be altered.\n",
                      msg, paste(changed, collapse = " and "))
        
        if (!is.null(matrixUse)) {
            ## Show stats for raw matrix
            ## Numeric conversion to prevent integer overflow on product
            rnr <- as.numeric( sum(populatedRows(matrixRaw)) )
            rnc <- as.numeric( sum(populatedCols(matrixRaw)) )
            rnz <- nnZero(matrixRaw)
            rpp <- smallNumberFormatter( rnz / (rnr * rnc) )
            ## Note the percent change for Use/Raw
            rcD <- doCol(sprintf("%.1f%%x%.1f%%", 100 * nr/rnr,
                                    100 * nc / rnc), "yellow")
            rzD <- smallNumberFormatter( nz / rnz )
            msg <- sprintf("%s  %s : %d x %d (%s), %d non-zero (%s; %s survive filters)\n",
                           msg, doCol("[Raw Matrix]", "red"), rnr, rnc, rcD,
                           rnz, doCol(rpp, "red"), doCol(rzD, "yellow"))
        }
        mdRow <- ifelse(is.null(matrixMD), 0, nrow(matrixMD))
        if (mdRow != 0 && !compact) {
            ## Sample metadata values
            mcols <- metadata_keys()
            msg   <- paste(c(msg, doCol(sprintf("  %d IDs have metadata assigned in up to %d keys, eg:\n", mdRow, length(mcols)), "blue"), collapse = ""))
            for (mcol in mcols) {
                ## Get a non-NA sample ID with this annotation
                notNA <- !is.na(matrixMD[[ mcol]])
                rowNum <- which(notNA)[1]
                ## OMG ITS SO PAINFUL JUST TO GET A VECTOR
                exid  <- as.character(matrixMD[ rowNum, "id", with=FALSE][[1]])
                exval <- matrixMD[ exid, mcol, with=FALSE][[1]]
                mline <- sprintf('    %s$metadata(key="%s", id="%s") # %s\n',
                                 whtName, doCol(mcol,"red"),
                                 doCol(exid, "cyan"),
                                 doCol(exval, "magenta"))
                msg <- paste(c(msg, mline), collapse ="")
            }
        }
        if (!compact) msg <- paste(msg, sprintf("%s$help() %s\n", whtName,
              doCol("# Get more help about this object", "yellow")),
                             collapse='', sep='')
        ## Left pad if requested
        if (CatMisc::is.def(pad)) msg <-
            paste(c(sprintf("%s%s", pad, unlist(base::strsplit(msg, "\n"))),
                    "\n"), collapse = "\n")
        ## Ending up with extraneous newlines at end
        sub("\n+$", "\n", msg)
    },

    help = function (color=NULL) {
        "\\preformatted{
Show compact help information about the object
      color - Should output be colorized. Default NULL, which checks
              $color()
}"
        if (is.null(color)) color <- useColor() # Use EventLogger setting
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        objName <- .self$.selfVarName("myMatrix")
        whtName <- doCol(objName, "white")
        comCol  <- "yellow"

        txt <- sprintf("
?AnnotatedMatrix   %s
%s                 %s
str(%s, max.lev=3) %s

%s
%s$map( inputIDs ) %s

%s
%s$autoFilter()                %s
%s$filterByScore(min,max)      %s
%s$filterByFactorLevel(levels) %s
%s$filterByMetadata(key,val)   %s
%s$rNames(rowNames)            %s
%s$cNames(colNames)            %s
%s$filterSummary()             %s
%s$appliedFilters()            %s
%s$reset()                     %s

%s
%s$metadata_keys()   %s         
%s$metadata(id, key) %s

",

doCol("# Built-in documentation on the class", comCol),
whtName, doCol("# Summary report of the object", comCol),
whtName, doCol("# Inspect the object structure", comCol),

doCol("### Primary operation",comCol),
whtName, doCol("# Pivot inputIDs from one matrix dimension to other", comCol),

doCol("### Filtering the matrix",comCol),
whtName, doCol("# Apply default filters (from file)", comCol),
whtName, doCol("# Remove cells failing min and/or max score limits", comCol),
whtName, doCol("# Keep (or discard) only certain factor levels", comCol),
whtName, doCol("# Discard rows/cols with matching metadata values", comCol),
whtName, doCol("# Restrict/reorder rownames", comCol),
whtName, doCol("# Restrict/reorder colnames", comCol),
whtName, doCol("# Summarize number of removed rows/cols", comCol),
whtName, doCol("# Text-representation of filters applied so far", comCol),
whtName, doCol("# Reset all filters", comCol),

doCol("### Metadata methods", comCol),
whtName, doCol("# Vector of all metadata keys (tagnames)", comCol),
whtName, doCol("# Query metadata for ids and/or keys", comCol)

)

        txt <- c(txt, doCol("### Other object methods\n", comCol))
        cmds <- list(
            showParameters = "Show currently set parameters",
            showLog        = "Show events and timing",
            nnZero         = "Number of non-zero cells",
            rCounts        = "Row counts = non-zero columns per each row",
            cCounts        = "Col counts = non-zero rows per each column",
            rNames         = "With no parameters, show current rownames",
            cNames         = "With no parameters, show current colnames",
            populatedRows  = "Logical vector of rows with non-zero data",
            populatedCols  = "Logical vector of columns with non-zero data",
            removeEmptyRows = "Prune matrix to remove empty rows",
            removeEmptyCols = "Prune matrix to remove empty cols",
            removeEmpty    = "Remove both empty rows and cols",
            'is.factor'    = "TRUE if matrix is treated as a factor",
            'as.gmt'       = "Dump current matrix to GMT text format",
            color          = "Toggle use of crayon (colorized output)",
            show           = "Pretty-print function, auto-invoked by object"
            )
        for (cmd in names(cmds)) {
            txt <- c(txt, sprintf("%s$%s() %s\n", whtName, cmd,
                                  doCol(paste("#",cmds[[cmd]]), comCol)))
        }
        txt <- c(txt, "\n", doCol("### Object fields\n", comCol))
        fields <- list(
            file       = "Path to the matrix source file",
            matrixRaw  = "Un-filtered sparse matrix as loaded from file",
            matrixUse  = "Sparse matrix after filtering (null if no filters)",
            lvlVal     = "Level names for factor matrices",
            filterLog  = "All filtered-out rows/columns, as a data.frame",
            matrixMD   = "Metadata data.table, both rows and columns",
            rowChanges = "Named vector of any row names that needed alteration",
            colChanges = "Named vector of any col names that needed alteration"
            )
        for (field in names(fields)) {
            txt <- c(txt, sprintf("%s$%s %s\n", whtName, field,
                                  doCol(paste("#",fields[[field]]), comCol)))
        }
        message(paste(txt, collapse='', sep=''))
        invisible(NULL)
    }
    
)
