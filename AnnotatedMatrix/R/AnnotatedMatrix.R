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
#' AnnotatedMatrix( help=TRUE ) # Show this help
#'
#' ## In general:
#' myAnnMat <- AnnotatedMatrix( file=NA, params=NA, autofilter=TRUE, ...)
#' ## Specific toy example matrix:
#' myAnnMat <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#'
#' myAnnMat$help()                 # High-level help
#' myAnnMat$ANYMETHOD( help=TRUE ) # Each method has detailed help available
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
#' @importFrom reshape2 dcast
#' @import Matrix
#' @importFrom CatMisc is.def is.something is.empty.field methodHelp
#' @import ParamSetI
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
                contains = c("ParamSetI")
                )

## Names for return values on filters
.filterNames <- c("CellsFiltered","RowsFiltered","ColsFiltered")

AnnotatedMatrix$methods(
    
    initialize = function(file=NA, params=NA, autofilter=TRUE,
                          help=FALSE, ... ) {
        "Create a new AnnotatedMatrix object; Invoke with AnnotatedMatrix(...)"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }
        callSuper( params=params, paramDefinitions="
Name        [character] Short Name assigned to the matrix
Description [character] Description for the matrix
ScoreDesc   [character] Describes what the matrix values (scores) represent
Source      [character] Primary source, presumably a URL
Authority   [character] The name of the authority responsible for the data

RowDim      [character] Name for the row dimension
ColDim      [character] Name for the column dimension
RowUrl      [character] Base URL for row names (%s placeholder for name)
ColUrl      [character] Base URL for column names (%s placeholder for name)

MinScore    [numeric] Minimum score recognized by $autoFilter()
MaxScore    [numeric] Maximum score recognized by $autoFilter()
MinRowCount [character] Minimum assignments per row, can be a percentage, recognized by $autoFilter()
MaxRowCount [character] Maximum assignments per row, can be a percentage, recognized by $autoFilter()
MinColCount [character] Minimum assignments per column, can be a percentage, recognized by $autoFilter()
MaxColCount [character] Maximum assignments per column, can be a percentage, recognized by $autoFilter()
KeepLevel   [character] List of preserved factor levels recognized by $autoFilter()
TossLevel   [character] List of discarded factor levels recognized by $autoFilter()
TossMeta    [character] Metadata value filter recognized by $autoFilter()
", ...)
        if (!CatMisc::is.def(file))
            err("AnnotatedMatrix must define 'file' when created", fatal = TRUE)
        file    <<- file
        fromRDS <<- FALSE
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
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        x <- c(0L, 0L, 0L)
        min <- param("MinScore")
        if (CatMisc::is.something(min)) {
            x <- x + filterByScore(min=min, reason=attr(min,"comment"))
        }
        max <- param("MaxScore")
        if (CatMisc::is.something(max)) {
            x <- x + filterByScore(max=max, reason=attr(max,"comment"))
        }
        for (kt in c("Keep", "Toss")) {
            isKeep <- kt == "Keep"
            if (is.factor()) {
                ## KeepLevel / TossLevel
                ktLvl <- param(sprintf("%sLevel", kt))
                if (CatMisc::is.something(ktLvl)) {
                    x <- x + filterByFactorLevel(ktLvl, keep=isKeep,
                                                 reason=attr(ktLvl,"comment"))
                }
            }
            meta <- param(sprintf("%sMeta", kt))
            if (CatMisc::is.something(meta)) {
                ## "Description LIKE {Deprecated}"
                ## "Type Not RegExp Pseudo.Gene"
                fbms <- CatMisc::parenRegExp('^(\\S+)\\s+(not\\s+)?(\\S+)\\s+(.+)', meta, unlist=FALSE)
                for (fbm in fbms) {
                    isNot  <- ifelse(fbm[2] == "", FALSE, TRUE)
                    shouldKeep <- kt == !isNot
                    v  <- fbm[4]
                    ## ToDo - smuggle in Row/Col at front of val
                    x <- x + filterByMetadata(key=fbm[1], val=v, type=fbm[3],
                                              keep=shouldKeep,
                                              reason=attr(meta,"comment"))
                }
            }
        }

        ## Count filters
        for (mm in c("Min", "Max")) {
            for (rc in c("Row", "Col")) {
                ## eg "MinColCount"
                filtCnt <- param(sprintf("%s%sCount", mm, rc))
                if (!CatMisc::is.something(filtCnt)) next
                r <- attr(filtCnt,"comment")
                if (mm == 'Min') {
                    x <- x + filterByCount(min=filtCnt, MARGIN=rc, reason=r)
                } else {
                    x <- x + filterByCount(max=filtCnt, MARGIN=rc, reason=r)
                }
            }
        }

        if (x[1] > 0 && recursive) x <- x + autoFilter( verbose=FALSE )
        if (x[1] > 0 && verbose) message(c("Automatic filters have masked",x[1],
              "cells,",x[2],"rows, and",x[3],"cols"), prefix="[-]",
              bgcolor='cyan', color='yellow')
        invisible(setNames(x, .filterNames))
    },

    .detailZeroedRowCol = function( obj, fail, metric, reason, help=FALSE) {
        "Internal method to tally rows and columns zeroed out during filtering"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        ## What rows / columns were affected by this operation? i=row, j=col
        ij   <- which( fail )
        ## Which triples are currently nonzero?
        nonZ <- obj@x != 0

        ## The steps below (find non-zero triples, find impacted
        ## triples, identifiy impacted with no survivors) are chosen
        ## because other approaches with vapply were very slow.
        
        ## Which rows have at least one nonzero entry?
        iOk  <- unique(obj@i[ nonZ ])
        ## Which rows were impacted by the current operation?
        iHit <- unique( obj@i[ ij ] )
        ## Remove the ones that still have an active row somewhere:
        iBye <- setdiff(iHit, iOk)
        iNum <- length(iBye)
        if (iNum > 0) {
            ## Some rows are now all zero because of the filter
            ## +1 -> Because Matrix Row/Col are zero-indexed!
           ids <- rownames(obj)[ iBye + 1 ]
            .filterDetails(id=ids, type="Row", metric=metric, reason=reason)
        }
        
        ## Which cols have at least one nonzero entry?
        jOk  <- unique(obj@j[ nonZ ])
        ## Which cols were impacted by the current operation?
        jHit <- unique( obj@j[ ij ] )
        ## Remove the ones that still have an active row somewhere:
        jBye <- setdiff(jHit, jOk)
        jNum <- length(jBye)
        if (jNum > 0) {
            ## Some cols are now all zero because of the filter
            ## +1 -> Because Matrix Row/Col are zero-indexed!
            ids <- colnames(obj)[ jBye + 1 ]
            .filterDetails(id=ids, type="Col", metric=metric, reason=reason)
        }
        c(iNum, jNum) # Return counts of zeroed-out rows,cols
    },

    filterByScore = function( min=NA, max=NA, filterEmpty=FALSE, reason=NA,
                             help=FALSE) {
        "Filter matrix cells by simple numeric tests"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        obj   <- matObj()
        rv    <- c(0L, 0L, 0L)
        ttxt  <- ""
        if (CatMisc::is.def(min)) {
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
                ttxt  <- paste("score <", min)
                ijz   <- .detailZeroedRowCol( obj, fail, ttxt, reason )
                rv    <- rv + c(numZ, ijz)
            }
            .addAppliedFilter("SCORE", min, reason, '<')
        }
        if (CatMisc::is.def(max)) {
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
                tt   <- paste("score >", max)
                ttxt <- if (ttxt == "") { tt } else { paste(ttxt,"or >",max) }
                ijz  <- .detailZeroedRowCol( obj, fail, tt, reason )
                rv   <- rv + c(numZ, ijz)
            }
            .addAppliedFilter("SCORE", max, reason, '>')
        }
        if (rv[1] != 0) {
            matrixUse <<- obj
            if (filterEmpty) {
                ## Strip empty rows and columns
                if (CatMisc::is.something(reason)) ttxt <- paste(ttxt, reason)
                removeEmpty(ttxt)
            }
        }
        invisible(setNames(rv, .filterNames))
    },

    filterByCount = function(MARGIN, min=NULL, max=NULL, relative=TRUE,
                             filterEmpty=FALSE, reason=NA, help=FALSE) {
        "Filter rows or columns based on number of non-zero connections"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        rv     <- c(0L, 0L, 0L)
        counts <- NULL
        reason <- reason[1]
        ## Which dimension is being filtered?
        if (grepl('(1|row)', MARGIN, ignore.case=TRUE)) {
            MARGIN <- 'Row'
            counts <- rCounts()
        } else if (grepl('(2|col)', MARGIN, ignore.case=TRUE)) {
            MARGIN <- 'Col'
            counts <- cCounts()
        } else {
            err("$filterByCount(MARGIN) should be one of 1,2,row,col",
                fatal=TRUE)
        }

        ## `counts` is generic at this point, it is either for "row"
        ## or "col" depending on `MARGIN`. I will be using "element"
        ## to generically refer to the 'things' in counts

        ## If relative is true, then we will consider '%' requests
        ## against the number of elements that are CURRENTLY
        ## populated. Otherwise we will consider all elements from the
        ## dimension, whether populated or not. Note that this can
        ## make a significant difference during recursive filtering
        denom <- if (relative) { sum(counts != 0) } else { length(counts) }
        obj   <- matObj()
        ttxt  <- ""
        
        if (!is.null(min)) {
            x <- normalizePercent(min[1], denom)
            ## Elements below the request
            fail <- counts < x & counts != 0
            numE <- sum(fail)
            if (numE > 0) {
                ## At least some elements have failed due to the filter
                ttxt <- sprintf("%s count < %s", MARGIN, min)
                inds <- which(fail) - 1 # indicies in Matrix are zero-indexed
                numZ <- sum( counts[ fail ] ) # Number of cells impacted
                ## Identify the actual cells being zeroed out:
                cells <- if (MARGIN == 'Row') {
                    is.element(obj@i, inds) & obj@x != 0
                } else {
                    is.element(obj@j, inds) & obj@x != 0                    
                }
                ## Safety check
                if (sum(cells) != numZ) err("Sanity count mismatch for filterByCount", prefix="[CODE ERROR]")
                obj@x[ cells ] <- 0
                ijz <- .detailZeroedRowCol( obj, cells, ttxt, reason )
                rv  <- rv + c(numZ, ijz)
            }
            .addAppliedFilter("COUNT", max, reason, paste(MARGIN, '>'))
        }

        if (!is.null(max)) {
            x <- normalizePercent(max[1], denom)
            ## Elements above the request
            fail <- counts > x & counts != 0
            numE <- sum(fail)
            if (numE > 0) {
                ## At least some elements have failed due to the filter
                tt    <- sprintf("%s count > %s", MARGIN, max)
                ttxt  <- if (ttxt == "") { tt } else { paste(ttxt,"or >",max) }
                inds  <- which(fail) - 1 # indicies in Matrix are zero-indexed!
                numZ  <- sum( counts[ fail ] ) # Number of cells impacted
                ## Identify the actual cells being zeroed out:
                cells <- if (MARGIN == 'Row') {
                    is.element(obj@i, inds) & obj@x != 0
                } else {
                    is.element(obj@j, inds) & obj@x != 0                    
                }
                ## Safety check
                if (sum(cells) != numZ) err("Sanity count mismatch for filterByCount", prefix="[CODE ERROR]")
                obj@x[ cells ] <- 0
                ijz <- .detailZeroedRowCol( obj, cells, tt, reason )
                rv  <- rv + c(numZ, ijz)
            }
            .addAppliedFilter("COUNT", max, reason, paste(MARGIN, '<'))
        }
        
        if (rv[1] != 0) {
            matrixUse <<- obj
            if (filterEmpty) {
                if (CatMisc::is.something(reason)) ttxt <- paste(ttxt, reason)
                removeEmpty(ttxt)
            }
        }
        invisible(setNames(rv, .filterNames))
    },

    filterByFactorLevel = function( x, keep=TRUE, ignore.case=TRUE,
                                   filterEmpty=FALSE, reason=NA, help=FALSE) {
        "Zero-out cells in the matrix that fail factor criteria"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        if (!is.factor()) {
            err("Can not filterByFactorLevel() - matrix is not a factor")
            return( invisible(setNames(NA,'FilterError') ) )
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
            return( invisible(setNames(NA,'FilterError') ) )
        }
        if (length(xVal) == 0) {
            err("No valid levels provided to filterByFactorLevel()")
            return( invisible(setNames(NA,'FilterError') ) )
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
        ## The reported operator should be for what is NOT kept:
        op <- if (k) {"!="} else {"=="}
        ## The filter text describes failing cells
        ttxt <- sprintf("Level %s %s", op,
                        paste(lnames, collapse=', '))
        .addAppliedFilter("LEVELS", lnames, reason, op)
        
        ## Working with the @x value vector for the sparse
        ## matrix. Determine which values are 'unwanted'
        fail <- is.element(obj@x, xVal)
        ## If our selection is being kept, we need to invert the mask
        ## AND not consider entries already zero:
        if (keep) fail <- !fail & obj@x != 0
        numZ <- sum( fail )
        if (numZ > 0) {
            ## Zero-out the failing cells
            obj@x[ fail ] <- 0
            ijz <- .detailZeroedRowCol( obj, fail, ttxt, reason )
            rv  <- rv + c(numZ, ijz)
            matrixUse <<- obj
            if (filterEmpty) {
                if (CatMisc::is.something(reason)) ttxt <- paste(ttxt, reason)
                removeEmpty(ttxt)
            }
        }
        invisible(setNames(rv, .filterNames))
    },

    filterByMetadata = function(key, val, MARGIN=NULL, keep=TRUE,type="like",
                                op='or', filterEmpty=FALSE, reason=NA,
                                help=FALSE) {
        "Zero-out rows or columns that match metadata criteria"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        ## Normalize the requested key
        key <- key[1]
        x   <- is.element(tolower(metadata_keys()), tolower(key))
        x   <- which(x)
        if (length(x) == 0) {
            err(paste("$filterByMetadata(key='", key,
                      "') does not match any keys"))
            return( invisible(setNames(NA,'FilterError') ) )
        }
        if (length(x) > 1) {
            ## Try case insensitive?
            x   <- is.element(metadata_keys(), key)
            x   <- which(x)
            if (length(x) != 1) {
                err(c("filterByMetadata(key =", key,
                      ") matches multiple columns"))
                return( invisible(setNames(NA,'FilterError') ) )
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
                err("$filterByMetadata(MARGIN) should be one of NULL,1,2,row,col", fatal=TRUE)
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
        if (keep) type <- c("Not", type)
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
        
        if (keep) fail = !fail

        rv     <- c(0L, 0L, 0L)
        metric <- paste(metric, collapse=' ')
        ## No matches
        if (sum(fail) == 0) return( invisible(setNames(rv, .filterNames)) )

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
                          sprintf("[%s %s%s%s]", margTxt, op,
                                  ifelse(keep, ' Keep', ''),
                                  ifelse(ic, ' ignore.case', '')))


        if (rv[1] != 0) {
            matrixUse <<- obj
            if (filterEmpty) {
                ## Strip empty rows and columns
                if (CatMisc::is.something(reason)) metric <- paste(metric,reason)
                removeEmpty(metric)
            }
        }
        invisible(setNames(rv, .filterNames))
    },

    nnZero = function(obj=NULL, help=FALSE, ...) {
        "Return the count of non-zero elements in the matrix"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        if (is.null(obj)) obj <- matObj(...)
        Matrix::nnzero( obj )
    },

    rCounts = function(obj=NULL, help=FALSE, ...) {
        "Return the counts of non-zero columns for each row"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        if (is.null(obj)) obj <- matObj(...)
        Matrix::rowSums(obj != 0)
    },

    populatedRows = function(obj=NULL, help=FALSE, ...) {
        "Determine which rows have at least one non-zero assignment"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        rCounts( obj=obj, ...) != 0
    },

    removeEmptyRows = function(reason=NA, help=FALSE) {
        "Remove all rows that consist only of zeroes"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        obj     <- matObj()
        isEmpty <- !populatedRows()
        ## populatedRows() returns a named logical vector. Take names
        ## from there
        toss    <- names(isEmpty)[ isEmpty ]
        if (length(toss) != 0) matrixUse <<- obj[ !isEmpty, , drop=FALSE]
        .addAppliedFilter( "REMOVE_EMPTY", "Rows", reason)
        invisible(toss)
    },

    cCounts = function(obj=NULL, help=FALSE, ...) {
        "Return the counts of non-zero rows for each column"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        if (is.null(obj)) obj <- matObj(...)
        Matrix::colSums(obj != 0)
    },

    populatedCols = function(obj=NULL, help=FALSE, ...) {
        "Determine which columns have at least one non-zero assignment"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        cCounts( obj=obj, ...) != 0
    },

    removeEmptyCols = function(reason=NA, help=FALSE) {
        "Remove all rows that consist only of zeroes"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        obj     <- matObj()
        isEmpty <- !populatedCols()
        ## populatedCols() returns a named logical vector. Take
        ## names from there
        toss    <- names(isEmpty)[ isEmpty ]
        if (length(toss) != 0) matrixUse <<- obj[ , !isEmpty, drop=FALSE]
        .addAppliedFilter( "REMOVE_EMPTY", "Cols", reason)
        invisible(toss)
    },

    removeEmpty = function(reason=NA, help=FALSE) {
        "Remove all rows and columns that consist only of zeroes"
        if (help) return(CatMisc::methodHelp(match.call(), class(.self),
                                             names(.refClassDef@contains)))
        toss <- removeEmptyRows(reason)
        toss <- c(toss, removeEmptyCols(reason))
        invisible(toss)
    },

    map = function(input=NULL, via=NULL, format="data.frame",
                   ignore.case=TRUE, keep.best=FALSE,
                   column.func=max,
                   collapse=NULL, collapse.name=NULL, collapse.token=',',
                   collapse.score=NULL, integer.factor=NULL,
                   add.metadata=TRUE, input.metadata=FALSE, warn=TRUE,
                   append.to=NULL, append.col=1L, recurse=0,
                   in.name="Input", out.name="Output", prefix.metadata=TRUE,
                   help=FALSE
                   ) {
        "Map (convert) names from one dimension of the matrix to the other"
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
            unk <- inpNms[ match(unk, inp) ] # Restore user's case
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

        ## recLevel, depth and bait are mainly relevant for recursive expansion
        recLevel <- 0L     # Current recursion level, 0 = user input
        depth    <- character() # Hash of minimum depth an ID was observed at
        ## Set of novel IDs being used for this recursion cycle. The
        ## IDs being queried are the names, the text to display
        ## (case-preserved) are the values.
        bait     <- setNames( inpNms[ match(ids, inp) ], ids )
        while (recLevel <= recurse && length(bait) > 0) {
            newFound      <- character() # Will hold newly-found IDs
            baitIds       <- names(bait)
            depth[ baitIds ] <- recLevel # Set recursion depth for current bait
            for (id in baitIds) {
                ## Cycle through each input id, building columns as we go:
                idNm <- bait[id]  # User's name (before to.lower)
                ## Select the indices for matrix row(s) matching the id:
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
                ## Note any new IDs for next recursion round. The
                ## vector will be cleaned up prior to starting the
                ## next cycle:
                if (recurse != 0L) newFound <- c(newFound, names(vec))

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
            recLevel <- recLevel + 1L
            if (recurse != 0L) {
                ## Find novel Outptut IDs to use for next round of
                ## recursion. Make found output IDs unique and
                ## normalize case:
                bait <- unique(newFound)
                bait <- if (ignore.case) {
                    setNames(bait, tolower(bait))
                } else {
                    setNames(bait, bait)
                }
                ## Only use IDs we have not already used:
                bait <- bait[ base::setdiff( names(bait), names(depth) ) ]
                ## Exclude any that are not in the input dimension:
                bait <- bait[ intersect(names(bait), matIds) ]
                ## All the Rube Goldberg naming is a mechanism to
                ## manage case-insensitive matching of IDs while still
                ## preserving 'original' case for reporting.
            }
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
        
        rv <- data.frame(Input  = inpCol[sl], Output = outCol[sl],
                         stringsAsFactors=FALSE)
        ## Customize column names if requested
        inName <- in.name[1]
        if (is.null(inName)) {
            ## Try to get name from dimensions
            d <- if (via == 'row') { 1 } else { 2 }
            inName <- names(dimnames(obj))[d]
            if (!CatMisc::is.something(inName)) inName <- "Input"
        }
        outName <- out.name[1]
        if (is.null(outName)) {
            ## Try to get name from dimensions
            d <- if (via == 'row') { 2 } else { 1 }
            outName <- names(dimnames(obj))[d]
            if (!CatMisc::is.something(outName)) outName <- "Output"
        }
        if (inName  != "Input")  names(rv)[1] <- inName
        if (outName != "Output") names(rv)[2] <- outName
        if (inName == outName) {
            err(c("Input and Output column headers are both", inName,
                  "- this will not end well, halting"))
            return(NA)
        }

        colHelp <- list()
        colHelp[[ inName ]]  <- "The input IDs provided to the $map() method"
        colHelp[[ outName ]] <- "The output IDs discovered from the input"

        depName <- NULL
        if (recurse != 0L) {
            ## Add a depth column
            depName <- "Depth"
            nms     <- rv[[ inName ]]
            if (ignore.case) nms <- tolower(nms)
            rv[[ depName ]] <- depth[ nms ]
            colHelp[[ outName ]] <- "The recursion depth the *INPUT* ID was first observed at"
        }
        
        
        if (!isFac || all(integer.factor)) {
            ## Add a numeric score column
            if (isFac) {
                rv$Score <- as.integer(scrCol[sl])
            } else {
                rv$Score <- scrCol[sl]
            }
            colHelp$Score <- "The recorded numeric value that connects input to output"
        }

        multOut <- base::duplicated(rv[[ outName ]])
        if (sum(multOut) == 0) {
            multOut <- character()
        } else {
            ## Some of the output terms come from multiple inputs
            inds    <- which(multOut)
            ## Note them so they can be reported to user:
            multOut <- unique( rv[[ outName ]][multOut] )
            if (colOut) {
                ## Request to collapse by input ID
                for (nm in multOut) {
                    ## Get all rows for this Output name:
                    nInds <- which(rv[[ outName ]] == nm)
                    ## We will keep the first row:
                    keep  <- nInds[1]
                    ## Aggregate the scores, stuff into the first index:
                    rv[keep, "Score"] <- valueCollapseFunction(scrCol[ nInds ])
                    ## ... and the input names:
                    rv[keep, inName]  <- nameCollapseFunction( inpCol[ nInds ])
                    if (recurse != 0L) {
                        ## Also collapse the recursion depth
                        rv[keep, depName] <-
                            nameCollapseFunction( rv[nInds, depName] )
                    }
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
            rv[[ outName ]]
        } else if (colIn || length(multIn) == 0) {
            ## Collapsing on Input, or Input (by chance or design) is unique:
            rv[[ inName ]]
        } else {
            strict.unique( rv[[ inName ]] )
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
            colHelp$Factor <- "The factor levels associated with the connection between input and output"
        }

        metacols <- NULL
        if (CatMisc::is.something(add.metadata)) {


### TODO: Need to manage cases where colIn is true, and Output is a
### token-separated string of IDs. That is, need to grab metadata in
### some fashion for these synthetic values. Should probably hold a
### temp list for each row of original values
            metacols <- character()
            metasrc  <- character()
            
            srcCols <- outName
            if (input.metadata) srcCols <- c(srcCols, inName)
            
            for (src in srcCols) {
                ## Request to include metadata columns
                ## Map our output IDs to the Metadata data.table :
                md  <- matrixMD[rv[[ src ]], setdiff(colnames(matrixMD), "id"),
                                with=FALSE]
                mdc <- colnames(md) # All available metadata columns
                ## If the param is not a character vector, take all
                if (!is.character(add.metadata)) add.metadata <- mdc
                for (col in add.metadata) {
                    rvCol <- col
                    if (prefix.metadata) rvCol <- paste(src, col)
                    nowNames <- colnames(rv)
                    if (is.element(rvCol, nowNames)) {
                        ## The metadata column will conflict with a column
                        ## already present. Make it unique
                        rvCol <- make.unique(c(nowNames, rvCol))[ length(nowNames) + 1 ]
                    }
                    colHelp[[ rvCol ]] <- paste("Metadata values for", col, "associated with", src)
                    metacols <- c(metacols, rvCol)
                    metasrc  <- c(metasrc, src)
                    if (is.element(col, mdc)) {
                        ## This is a known metadata column. Include it,
                        ## unless it's all NAs
                        if (!all(is.na(md[[ col ]]))) rv[[rvCol]] <- md[[ col ]]
                    } else {
                        ## Column does not exist, put a blank text column
                        ## in the output
                        rv[[rvCol]] <- character(rcnt)
                    }
                }
            }
            attr(metacols, "Source") <- metasrc
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
            rv <- setNames(rv[[ outName ]], rv[[ inName ]])
            ## This can result in non-unique names. Do we want to use
            ## strict.unique() to uniquify names?
            class(rv) <- c("mapResult", class(rv))
        } else if (grepl('matrix', format[1], ignore.case=TRUE)) {
            ## Matrix format
            ## Find all valid edges (both Input and Output defined):
            ok <- !is.na(rv[[ inName ]]) & !is.na(rv[[ outName ]])
            sc <- rv$Score
            if (is.null(sc)) {
                err("Score column is NULL, all weights set to 1 in matrix")
                sc <- rep(1L, sum(ok))
            } else {
                ok <- ok & !is.na(sc)
                sc <- sc[ ok ]
            }
            
            rn <- unique(rv[ ok, inName ])
            cn <- unique(rv[ ok, outName ])

            i <- match(rv[ ok, inName ], rn)
            j <- match(rv[ ok, outName ], cn)
            dn <- list()
            dn[[ inName  ]] <- rn
            dn[[ outName ]] <- cn

### Not working when consumed by qgraph(). Weird resulting graph,
### initially looks reasonable, but has incorrect connections.  I
### think qgraph may not handle sparse Matrices? I think the
### underlying matrix is accurate (TODO: need more tests)
            
            rv <- as(sparseMatrix( i=i, j=j, x=sc, dimnames=dn), "dgTMatrix")
            
        } else if (grepl('dynamic', format[1], ignore.case=TRUE)) {
            ## dynamictable HTML output
            if (require("dynamictable", quietly=TRUE)) {
                ord   <- c('row', 'col')
                if (via == 'col') ord <- rev(ord)
                urlParams <- unlist(lapply(ord, function(x) {
                    param(sprintf("%surl", x))
                }))
                
                opts  <- list( )
                if (isFac) {
                    opts$Factor = list(factor=TRUE)
                    opts$Score  = list(byFactor="Factor")
                }
                opts[[ inName ]]  <- list( url = urlParams[1] )
                opts[[ outName ]] <- list( url = urlParams[2] )
                opts[["*"]] <- list(truncate=80)
                
                rv <- dynamictable( rv, options=opts, auto.title=FALSE )
            } else {
                err("The dynamictable package is not installed; Please run install.packages('dynamictable')")
            }
        } else if (grepl('canvas', format[1], ignore.case=TRUE)) {
            ## CanvasXpress DHTML network format
            if (require("canvasXpress", quietly=TRUE)) {
                ## Exclude rows that don't have both input and output:
                notEdge  <- is.na(rv[[inName]]) | is.na(rv[[outName]])
                inNode   <- rv[ !notEdge, inName ]
                outNode  <- rv[ !notEdge, outName ]
                ## Get all unique nodes from both input and output
                nNames   <- unique(c(inNode, outNode))
                ## Capture requested metadata associated with the nodes:
                nodes    <- metadata( nNames , key=add.metadata, drop=FALSE )
                ## Build basic edge structure:
                edges    <- data.frame(id1=inNode, id2=outNode,
                                       stringsAsFactors=FALSE)
                outCols  <- colnames(rv)
                portCols <- c("Depth", "Score", "Factor")
                for (pc in portCols) {
                    if (is.element(pc, outCols))
                        edges[[ pc ]] <- rv[ !notEdge, pc]
                }
                title <- sprintf("Map from %s", if (length(inp) <= 3) {
                    paste(inp, collapse=' ')
                } else {
                    paste(length(inp), "nodes")
                })
                rv <- canvasXpress(
                    nodeData=nodes,
                    edgeData=edges,
                    edgeWidth=2,
                    graphType="Network",
                    nodeFontColor="rgb(29,34,43)",
                    nodeSize=30,
                    showAnimation=FALSE,
                    directed=TRUE,
                    ## networkLayoutType="radial",
                    ## networkRoot=root,
                    title=title )
            } else {
                err("The canvasXpress package is not installed; Please run install.packages('canvasXpress')")
            }
        } else if (grepl('graph', format[1], ignore.case=TRUE)) {
            ## graph format
            if (require("qgraph", quietly=TRUE)) {
                nodes <- unique(c(rv[[inName]], rv[[outName]]))
                pairs <- matrix(c(match(rv[[inName]], nodes),
                                  match(rv[[outName]], nodes)), ncol=2)
                ## I was pretty sure there was a way to linearize a matrix
                ## by row, but I can't remember it. So transpose the
                ## In/Out matrix and vectorize:
                edges <- as.vector(t(pairs))
                g <- graph(edges, directed=TRUE)
                vertex_attr(g, "name") <- nodes
                if (!is.null(rv$Score))  edge_attr(g, "weight") <- rv$Score
                if (!is.null(rv$Factor)) edge_attr(g, "name")   <- rv$Factor
                rv <- g
            } else {
                err("The igraph package is not installed; Please run install.packages('igraph')")
            }
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
        attr(rv, "Columns")    <- colHelp
        attr(rv, "Metadata")   <- metacols
           
        if (warn) {
            ## Alert user to any non-unique relationships:
            msg <- character()
            cm  <- function(w,x,cl) colorize(paste(w,":",length(x)), cl)
            if (length(unk) != 0)     msg <- c(msg,cm("Unknown",unk,"magenta"))
            if (length(unMap) != 0)   msg <- c(msg,cm("Unmapped",unMap,"red"))
            if (length(dupIn) != 0)   msg <- c(msg,cm("Dup.In",dupIn,"blue"))
            if (length(dupMat) != 0)  msg <- c(msg,cm("Dup.Mat",dupMat,"cyan"))
            if (length(multIn) != 0)  msg <- c(msg,cm("Mult.In",multIn,"yellow"))
            if (length(multOut) != 0) msg <- c(msg,cm("Mult.Out",multOut,"green"))
            if (length(msg) !=0) base::message("Non-unique events: ",
                                               paste(msg, collapse=', '))
        }

        rv
    },

    .autoLevel = function (vals, sep=',', decreasing=TRUE, help=FALSE) {
        "Generate new factor levels given a set of existing ones"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
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

    .filterDetails = function (id=NA, type=NA,  metric=NA, reason=NA,
                               help=FALSE) {
        "Utility method to extend filter log"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
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

    .addAppliedFilter = function (key, val, com=NULL, pre=NULL, pro=NULL,
                                  help=FALSE) {
        " Internal method, stores filter text in $setFilters"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        txt <- c(toupper(key))
        if (CatMisc::is.something(pre)) txt <- c(txt, pre)
        txt <- c(txt, paste(val, collapse=sfSep))
        if (CatMisc::is.something(pro)) txt <- c(txt, pro)
        if (CatMisc::is.something(com)) txt <- c(txt, "##", com)
        txt <- paste(txt, collapse=' ')
        if (!is.element(txt, setFilters)) setFilters <<- c(setFilters, txt)
        txt
    },

    appliedFilters = function (new=NULL, help=FALSE) {
        "Get applied filters as a vector of readable/parsable strings"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!is.null(new)) {
            err("
ToDo: STILL WORKING ON ROUND-TRIP PARSING FILTER TEXT
")
        }
        setFilters
    },

    filterSummary = function (reason=TRUE, help=FALSE) {
        "Human-readable overview of filters and counts of rows/cols removed"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )

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

    is.factor = function (help=FALSE) {
        "TRUE if the matrix is a pseudo-factor, otherwise FALSE"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                              names(.refClassDef@contains)) )
        CatMisc::is.something(lvlVal)
    },

    levels = function( asFactor=FALSE, help=FALSE ) {
        "Returns factor levels, if appropriate. If not, returns NULL"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (is.factor()) {
            if (asFactor) { factor(lvlVal) } else { lvlVal }
        } else {
            NULL
        }
    },

    as.gmt = function( obj=NULL, transpose=FALSE, file=NULL, help=FALSE, ... ) {
        "Convert the active matrix into GMT text representation"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
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

    .readMatrix = function ( format = "", help=FALSE, ... ) {
        "Internal method to read matrix data from a variety of formats"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
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

    metadata_keys = function ( help=FALSE ) {
        "Get all metadata keys as a character vector"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        ## Exclude the id column
        base::setdiff(names(matrixMD), "id")
    },

    metadata = function ( id=NULL, key=NULL, na.rm=TRUE, drop=TRUE,
                          verbose=TRUE, help=FALSE) {
        "Recover metadata for specific IDs and/or columns"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                              names(.refClassDef@contains)) )
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
                naCol <- sapply(useKey, function(x) is.na( rv[[ x ]] ),
                                simplify='matrix')
                ## And now find rows that are all() NAs:
                naRow <- if (is.matrix(naCol)) {
                    ## sapply honored our request to stay a matrix
                    apply(naCol, 1, all)
                } else if (nrow(rv) <= 1) {
                    ## rv was 0-1 rows, use it as the mask
                    all(is.na(naCol))
                } else {
                    ## So we must have had a single column, the vector
                    ## is our mask
                    naCol
                }
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
                naCol1 <- sapply(cols, function(x) is.na( rv[[ x ]] ),
                                simplify='matrix' )
                ## And then columns that are all NAs:
                naCol2 <- if (is.matrix(naCol1)) {
                    ## sapply honored our request to stay a matrix
                    apply(naCol1, 2, all)
                } else if (nrow(rv) <= 1) {
                    ## rv was 0-1 rows, use it as the mask
                    naCol1
                } else {
                    ## So we must have had a single column
                    all(is.na(naCol1))
                }
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

    .fieldDescriptions = function ( update=TRUE, help=FALSE ) {
        "Update object fields to include attributes with brief descriptions"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                              names(.refClassDef@contains)) )
        fields <- list(
            file       = "Path to the matrix source file",
            fromRDS    = "Logical flag indicating file was read from cached RDS",
            matrixRaw  = "Un-filtered sparse matrix as loaded from file",
            matrixUse  = "Sparse matrix after filtering (null if no filters)",
            matrixMD   = "Metadata data.table, both rows and columns",
            filterLog  = "All filtered-out rows/columns, as a data.frame",
            setFilters = "Machine- and human-parsable list of applied filters",
            lvlVal     = "Level names for factor matrices",
            rowChanges = "Named vector of any row names that needed alteration",
            colChanges = "Named vector of any col names that needed alteration"
        )
        if (update) {
            hfmt <- " help('%s', 'AnnotatedMatrix') # More information on field "
            for (fld in names(fields)) {
                if (is.null(.self[[fld]])) next # Can't attribute NULL
                ## The [[ accessor seems to work for fields?
                attr(.self[[fld]], "Description") <- fields[[fld]]
                attr(.self[[fld]], "Help") <- sprintf(hfmt,fld)
            }
        }
        fields
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
        "Display high-level help about all object methods"
        if (is.null(color)) color <- useColor() # Use EventLogger setting
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        objName <- .self$.selfVarName("myMatrix")
        whtName <- doCol(objName, "white")
        comCol  <- "yellow"

        ## Figure out all available methods:
        allMeth <- AnnotatedMatrix$methods()
        ## Subtract out the generic ones:
        allMeth <- setdiff(allMeth, c("callSuper", "copy", "export", "field", "getClass", "getRefClass", "import", "initFields", ".objectPackage", ".objectParent", "show", "trace", "untrace", "usingMethods"))
        ## Subtract out some superclasses
        allMeth <- allMeth[ !grepl('#', allMeth) ]

        ## Organization of primary methods:
        sections <- list(
            "Primary Operation" = c("map", "as.gmt"),
            "Filtering the Matrix" = c("filterByScore", "filterByCount", "filterByFactorLevel","filterByMetadata", "rNames", "cNames", "removeEmptyRows", "removeEmptyCols", "removeEmpty", "reset", "autoFilter", "filterSummary", "appliedFilters"),
            "Matrix Information" = c("rCounts", "cCounts","populatedRows", "populatedCols", "nnZero", "is.factor", "levels"),
            "Metadata" = c("metadata", "metadata_keys"),
            "Parameter Management" = c("param", "showParameters", "allParams", "defineParameters","paramClass", "paramDefinition", "paramName", "setParamList", "hasParam" ),
            "Event Logging" = c("message","showLog","verbose","useColor" ),
            "SKIP" = c("actionMessage","dateMessage","debugMessage","err"),
            "Internal Methods" = allMeth[ grepl('^\\.', allMeth) ])
        ## Everything else:
        sections[["Other Methods"]] <- setdiff(allMeth, unname(unlist(sections)))

        ## Basic header:
        txt <- sprintf("
%s
?AnnotatedMatrix %s
%s                 %s
%s$showLog()       %s
str(%s, max.lev=3) %s
",
doCol("###
### Annotated Matrix Help - call the below commands for more details
###","magenta"),
doCol("# Built-in documentation on the class", comCol),
whtName, doCol("# Summary report of the object", comCol),
whtName, doCol("# Show logged events with timing", comCol),
whtName, doCol("# Inspect the object structure", comCol))

        noHelp <- c()
        ## Add snippets for each method, broken down by section
        for (sec in names(sections)) {
            if (sec == "SKIP") next
            txt <- c(txt, doCol(paste("\n############\n###", sec, "\n"),comCol))
            meths <- sections[[ sec ]]
            for (meth in meths) {
                ## Going to see if we can extract the ROxygen
                ## description string from the method
                code <- capture.output(AnnotatedMatrix$methods(meth))
                ## Should not happen, but be safe:
                if (is.null(code)) next
                isHelped <- FALSE
                com <- NA
                for (line in code) {
                    ## See if 'help=FALSE' is set - indicates I have
                    ## tied it into my internalized help framework:
                    if (grepl('help\\s*=\\s*FALSE', line))  isHelped <- TRUE
                    cm <- CatMisc::parenRegExp('^\\s+"(.+)"\\s*$', line)
                    if (!is.na(cm[1])) {
                        ## Found a single line quoted string, presume
                        ## it is description and stop scanning code:
                        com <- cm[1]
                        break
                    }
                }
                if (!isHelped) {
                    noHelp <- c(noHelp, meth)
                    next
                }
                txt <- c(txt, sprintf("%s$%s( help=TRUE )", whtName, meth))
                if (!is.na(com)) txt <-
                     c(txt, doCol(paste("\n    #", com), comCol))
                txt <- c(txt, "\n")
            }
        }
        if (length(noHelp) > 0) txt <- c(txt,
              doCol("\n### Methods lacking help\n", comCol),
              sprintf("# %s\n", strwrap(paste(noHelp, collapse=' '))))

        txt <- c(txt, doCol("\n### Object fields\n", comCol))
        fields <- .fieldDescriptions( )
        for (field in names(fields)) {
            txt <- c(txt, sprintf("str(%s$%s) %s\n", whtName, field,
                                  doCol(paste("#",fields[[field]]), comCol)))
        }
        base::message(paste(txt, collapse='', sep=''))
        invisible(NULL)
    }
)
