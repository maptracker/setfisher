#' Annotated Matrix
#'
#' Annotated sparse matrix for capturing query lists, identifier
#' mappings and ontology lookups
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
#' @field rowChanges Named character vector of any row names that
#'     needed changing. Values are the original name, names are the
#'     names after processing with make.names() (if valid = TRUE) or
#'     make.unique() if (valid = FALSE)
#' @field colChanges As per rowChanges, but for column names
#'
#' @importFrom methods new setRefClass
#' @importFrom utils read.table
#' @importFrom data.table data.table rbindlist as.data.table fread
#'     setkey setkeyv fsetdiff
#' @import Matrix
#' @importFrom CatMisc is.def is.something
#' @import ParamSetI
#' @importClassesFrom EventLogger EventLogger
#'
#' @examples
#'
#' ## Load a toy symbol-to-gene mapping matrix and use it to convert
#' ## some genes to Entrez Gene IDs
#' demo("geneSymbolMapping", package="AnnotatedMatrix", ask=FALSE)
#'
#' ## Work with different file formats
#' demo("fileFormats", package="AnnotatedMatrix", ask=FALSE)
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
                    lvlVal     = "character",

                    ## If row or column names need to be remapped:
                    rowChanges = "character",
                    colChanges = "character"
                    ),
                contains = c("EventLogger", "ParamSetI")
                )


AnnotatedMatrix$methods(
    
    initialize = function(file=NA, params=NA, ... ) {
        "\\preformatted{
Create a new object using AnnotatedMatrix():
       file - Required, a path to the file that stores the matrix. See
              .readMatrix() for details on supported formats
     params - Optional list of key/value pairs that will be passed to
              setParamList()
        ... - dots will also be passed on to setParamList(), as well as to
              .readMatrix()
}"
        callSuper(...)
        if (!CatMisc::is.def(file))
            err("AnnotatedMatrix must define 'file' when created", fatal = TRUE)
        file    <<- file
        fromRDS <<- FALSE
        defineParameters("
Name        [character] Optional name assigned to the matrix
Description [character] Optional description for the matrix
RowDim      [character] Optional name for the row dimension
ColDim      [character] Optional name for the column dimension
RowUrl [character] Optional base URL for row names (%s placeholder for name)
ColUrl [character] Optional base URL for column names (%s placeholder for name)
")
        setParamList( params=params, ... )
        .readMatrix( ... )
        reset()
    },
    
    matrix = function( raw=FALSE, transpose=FALSE) {
        "\\preformatted{
Retrieves the underlying Matrix for this object. Parameters:
      raw - Default FALSE, in which case the filtered Matrix (held in field
            'matrixUse') will be returned, if it is available. If not available,
            or if raw is TRUE, then the raw (as loaded from file) Matrix
            will be returned.
 transpose - Default FALSE. If true, the matrix will be transposed and names
            will be transfered to the appropriate dimensions
}"
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

    reset = function( asFactor=FALSE ) {
        "Reset any filters that were applied - the 'used' matrix will be the original 'raw' one"
        matrixUse <<- NULL
        filterLog <<- data.table(id = character(), key = "id")
        invisible(NA)
    },

    filterByScore = function( min=NA, max=NA, filterEmpty=FALSE, reason=NA ) {
        "\\preformatted{
Apply filters to the current matrix to zero-out cells failing thresholds.
        min - Minimum allowed cell value. Cells below this will be set to zero
        max - Maximum allowed cell value. Cells above it will be set to zero
 filterEmpty - Default FALSE; If true, then the matrix will be 'shrunk' to
              remove rows and columns that are only zeros
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        obj <- matrix()
        rv  <- 0
        if (filterEmpty) removeEmpty("Empty rows and cols before score filter")
        if (is.something(min)) {
            ## Zero out entries that fall below min
            inds <- which( obj != 0 & obj < min, arr.ind = TRUE)
            numZ <- nrow(inds)
            if (numZ > 0) {
                ## At least some cells were zeroed out
                obj[ inds ] <- 0
                rv      <- rv + numZ
                testTxt <- paste("x <", min)
                numTxt  <- paste(numZ, "Cells")
                filterDetails(id=testTxt, type="Val", metric=numTxt,
                              reason=reason)
                if (filterEmpty) {
                    ## Strip empty rows and columns
                    matrixUse <<- obj
                    if (!is.na(reason[1])) testTxt <- paste(testTxt, reason[1])
                    removeEmpty(testTxt)
                    obj <- matrix()
                }
            }
        }
        if (is.something(max)) {
            ## Zero out entries that are above max
            inds <- which( obj != 0 & obj > max, arr.ind = TRUE)
            numZ <- nrow(inds)
            if (numZ > 0) {
                ## At least some cells were zeroed out
                obj[ inds ] <- 0
                rv      <- rv + numZ
                testTxt <- paste("x >", max)
                numTxt  <- paste(numZ, "Cells")
                filterDetails(id=testTxt, type="Val", metric=numTxt,
                              reason=reason)
                if (filterEmpty) {
                    ## Strip empty rows and columns
                    matrixUse <<- obj
                    if (!is.na(reason[1])) testTxt <- paste(testTxt, reason[1])
                    removeEmpty(testTxt)
                    obj <- matrix()
                }
            }
        }
        if (rv != 0 && !filterEmpty) matrixUse <<- obj
        invisible(rv)
    },

    populatedRows = function(obj=NULL, ...) {
        "\\preformatted{
Return a logical vector indicating which rows have at least one non-zero cell
        obj - Default NULL, which will recover the matrix from matrix(),
              passing ... as well (so you can select the raw matrix if desired)
}"
        if (is.null(obj)) obj <- matrix(...)
        Matrix::rowSums(obj != 0) != 0
    },

    removeEmptyRows = function(reason=NA) {
        "\\preformatted{
Remove all empty rows (those that only contain zeros). Invisibly returns a
vector of removed IDs.
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        obj     <- matrix()
        isEmpty <- !populatedRows()
        toss    <- names(obj)[ isEmpty ]
        if (length(toss) != 0) {
            ## Some columns have been removed
            filterDetails(id=toss, type="Row", metric="AllZero", reason=reason)
            matrixUse <<- obj[ !isEmpty, ]
        }
        invisible(toss)
    },

    populatedColumns = function(obj=NULL, ...) {
        "\\preformatted{
Return a logical vector indicating which columns have at least one non-zero cell
        obj - Default NULL, which will recover the matrix from matrix(),
              passing ... as well (so you can select the raw matrix if desired)
}"
        if (is.null(obj)) obj <- matrix(...)
        Matrix::colSums(obj != 0) != 0
    },

    removeEmptyColumns = function(reason=NA) {
        "\\preformatted{
Remove all empty columns (those that only contain zeros). Invisibly returns a
vector of removed IDs.
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        obj     <- matrix()
        isEmpty <- !populatedColumns()
        toss    <- names(obj)[ isEmpty ]
        if (length(toss) != 0) {
            ## Some columns have been removed
            filterDetails(id=toss, type="Col", metric="AllZero", reason=reason)
            matrixUse <<- obj[ , !isEmpty]
        }
        invisible(toss)
    },

    removeEmpty = function(reason=NA) {
        "\\preformatted{
Remove all empty rows and columns. Invisibly returns a vector of removed IDs.
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        toss <- removeEmptyRows(reason)
        toss <- c(toss, removeEmptyColumns(reason))
        invisible(toss)
    },

    map = function(input, via=NULL, ignore.case=TRUE, column.func=max,
                   keep.best=FALSE,
                   collapse=NULL, collapse.name=NULL, collapse.token=',',
                   collapse.func=mean,
                   collapse.factor=NULL, integer.factor=FALSE,
                   add.metadata=TRUE, warn=TRUE
                   ) {
        "\\preformatted{
Provide a list of IDs, and map/pivot it from one dimension of the matrix to the
other, following 'connections' defined by non-zero cells. Returns a data.frame
with Input and Output columns, plus Score and/or Factor columns.
      input - Required, a vector of IDs
        via - Specify if the input matches the 'rows' or 'columns' of the
              matrix. If NULL (default) then your input will be compared to the
              row and column names, and the one with the most matches will be
              chosen (defaulting to 'row' in the event of equal matches)
 ignore.case - Default TRUE, which ignores the capitilazation of IDs
 column.func - Default max. If ignore.case is true, it is possible that an
              input ID can match multiple matrix IDs. In this case, multiple
              matching rows will be returned for one ID. column.func is
              applied to reduce this to a single row.
  keep.best - Default FALSE. If TRUE, then only the top-scored cell(s) will
              be kept
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
 collapse.func - A function, defaulting to mean, used to collapse scores to
              a single value. Can be any other function, including user-supplied
 collapse.factor - A function  that will be used after collapse.func if the
              matrix is a factor (levels have been set). The default is NULL,
              which will result in $.autoLevel() being used to generate new
              'hybrid' factors as needed, but can be a user-supplied function.
              The function should presumably generate an integer value that
              will correspond to a level.
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
}"
        if (!is.vector(input) || !is.character(input[1])) {
            err("mapToCol() must be provided with a character vector as input")
            return(NA)
        }
        obj   <- matrix()
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
        dupIn <- duplicated(inp)
        if (sum(dupIn) == 0) {
            dupIn <- character()
        } else {
            ## Some of the input is duplicated
            dupIn <- unique( inp[dupIn] )
            dupIn <- vapply(dupIn, function(x) {
                paste(names(rn)[which(rn == x)], collapse=collapse.token)
            }, "")
            inp  <- inp[ !duplicated(inp) ]
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
            obj <- matrix(transpose=TRUE) # transpose for simplicity below
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
        valueCollapseFunction <- if (isFac) {
            ## Factor-based functions
            if (is.null(collapse.factor)) {
                ## Default method to handle factors by making
                ## new 'hybrid' ones
                function (vec) .autoLevel( vec, sep=collapse.token )
            } else {
                ## User-supplied factor handling function
                collapse.factor
            }
        } else {
            ## Basic method to generate a single numeric value -
            ## default is mean()
            collapse.func
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
                rows  <- base::matrix(0L, 1, 1, dimnames=list(r=NA, c=NA))
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
                    colVal <- valueCollapseFunction(vec)
                    nmVal  <- nameCollapseFunction( names(vec) )
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

        multOut <- duplicated(rv$Output)
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
        
        if (isFac && !any(integer.factor)) {
            ## Add factor text column
            inds <- as.integer(scrCol[sl])
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

        if (is.something(add.metadata)) {


### TODO: Need to manage cases where colIn is true, and Output is a
### token-separated string of IDs. Should probably hold a temp list
### for each row of original values
            
            ## Request to include metadata columns
            ## Map our output IDs to the Metadata data.table :
            md  <- matrixMD[.(rv$Output), -"id"]
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

        ## Attach some summary attributes:
        attr(rv, "Unmapped")   <- unMap
        attr(rv, "Unknown")    <- unk
        attr(rv, "Via")        <- via
        attr(rv, "Dup.In")     <- dupIn
        attr(rv, "Dup.Mat")    <- dupMat
        attr(rv, "Mult.In")    <- multIn
        attr(rv, "Mult.Out")   <- multOut
        if (isFac) attr(rv, "Levels") <- lvlVal
        ## Also attach a brief explanation of each attribute:
        attr(rv, "Notes")      <-
            list(Via = "Whether your input was matched to rows or columns of the matrix",
                 'Dup.In' = "IDs that were present twice or more in input (possibly after case removal)",
                 'Dup.Mat' = "IDs that were present twice or more in matrix (after case removal)",
                 'Mult.In' = "Input IDs that generated multiple output values",
                 'Mult.Out' = "Output IDs that were generated from multiple inputs",
                 Unmapped = "Input ID is also in the matrix, but does not have a target with non-zero score",
                 Unknown = "Input ID could not be matched to any in the matrix")
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

    filterDetails = function ( id=NA, type=NA,  metric=NA, reason=NA) {
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
        row <- data.table::data.table(id=id, metric=metric, type=type,
                                      reason=reason, key="id")
        ## Do not add any entries that are already recorded as
        ## filtered (only note the first exclusion)
        newR <- setdiff(row[["id"]], filterLog[["id"]])
        filterLog <<- data.table::rbindlist(list(filterLog, row[newR]),
                                            fill = TRUE)
        data.table::setkeyv(filterLog, "id")
        filterLog
    },

    is.factor = function () {
        "TRUE if the matrix is a pseudo-factor (levels have been defined), otherwise FALSE"
        CatMisc::is.something(lvlVal)
    },

    levels = function( asFactor=FALSE ) {
        "\\preformatted{Returns factor levels, if appropriate. If not, returns NULL
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
        obj - Default NULL, which will recover the matrix from matrix(),
              passing ... as well. Alternatively a Matrix can be provided.
  transpose - Default FALSE, which will presume that the rows are sets. If TRUE,
              then columns will be taken as sets
       file - Default NULL, if defined then the output will be written to that
              path
}"
        ## http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
        if (is.null(obj)) obj <- matrix(transpose=transpose, ...)
        ## Only keep rows (sets) with at least one object:
        hasData  <- populatedRows(obj)
        obj      <- obj[hasData, ]
        setNames <- rownames(obj) # Rows are sets
        memNames <- colnames(obj) # Columns are the potential members
        descr    <- matrixMD[.(setNames), "Description"]
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

    ## Metadata can be held either in the comments of the MatrixMarket
    ## file, or as a separate file with a related name. Using a
    ## separate file is faster (roughly two fold)

    ## Separate .mtx and metadata files (~3.8 sec)

    ## 454.137 ms | Reading Matrix file
    ##              MAP_HG-U219_to_EntrezGene.mtx
    ##    2.251 s | Parsing Row names
    ## 978.887 ms | Parsing Col names
    ## 101.468 ms | Reading metadata file -
    ##              MAP_HG-U219_to_EntrezGene.mtx-Metadata.tsv
    ##            | Finished - 48867 x 20823, 0.00498% populated

    ## Unified .mtx file with embdedded metadata (~8.5 sec, 31Mb file)

    ##    1.049 s | Reading Matrix file
    ##              MAP_HG-U219_to_EntrezGene.mtx
    ##    5.251 s | Parsing Row names
    ##    2.286 s | Parsing Col names
    ##            | Finished - 48867 x 20823, 0.00498% populated

    ## Serializing timings (3.3Mb file)
    ##    1.569 s | Serializing matrix to file HG-U219LL-GO.mtx.rds
    ##    1.233 s | Reading serialized matrix HG-U219LL-GO.mtx.rds
 

    .readMatrix = function ( format = "", ... ) {
        if (!CatMisc::is.def(file)) err("AnnotatedMatrix objects must define 'file' when created", fatal = TRUE)
        objFile <- paste(file,'rds', sep = '.')
        if (file.exists(objFile) &&
            (!file.exists(file) || file.mtime(objFile) >= file.mtime(file))) {
            ## A serialized version of the data already exists, and
            ## the original file is either missing, or it is older
            ## than the serialized one (ie, the raw file has not been
            ## updated since the serialization occured)
            dateMessage(paste("Reading serialized matrix",
                              colorize(objFile,"white")), prefix = "  ")
            rv <- readRDS(objFile)
            fromRDS <<- TRUE
        } else {
            if (!file.exists(file)) err(c("Can not make AnnotatedMatrix",
                                          "File does not exist : ", file),
                                        fatal = TRUE)

            if (grepl('(mtx|matrixmarket)', format, ignore.case = TRUE) ||
                grepl('\\.mtx', file, ignore.case = TRUE)) {
                rv <- .readMatrixMTX( ... )
            } else if (grepl('(txt|text)', format, ignore.case = TRUE) ||
                       grepl('\\.(txt|text|list)', file, ignore.case = TRUE)) {
                rv <- .readMatrixTXT( ... )
            } else if (grepl('(lol)', format, ignore.case = TRUE) ||
                       grepl('\\.(inp)', file, ignore.case = TRUE)) {
                rv <- .readMatrixLOL( ... )
            } else if (grepl('(gmt)', format, ignore.case = TRUE) ||
                       grepl('\\.(gmt)', file, ignore.case = TRUE)) {
                rv <- .readMatrixGMT( ... )
            } else {
                err(c("Can not make AnnotatedMatrix - unrecognized file type: ",
                      file), fatal = TRUE)
            }
            ## Stub DT if no metdata was defined:
            if (!CatMisc::is.something(rv$metadata))
                rv$metadata <- data.table::data.table( id = character(), key = "id" )
            dateMessage(paste("Serializing matrix to file",
                              colorize(objFile,"white")), prefix = "  ")
            saveRDS(rv, objFile)
        }
        matrixRaw  <<- rv$matrix
        matrixMD   <<- rv$metadata
        if (!is.null(rv$rowChanges)) rowChanges <<- rv$rowChanges
        if (!is.null(rv$colChanges)) colChanges <<- rv$colChanges
        if (CatMisc::is.def(rv$levels)) lvlVal <<- rv$levels
        ## Set default parameters, without clobbering any already set
        if (CatMisc::is.def(rv$params)) setParamList(rv$params, clobber = FALSE)

        ## Numeric conversion to prevent integer overflow on product
        cnum       <- as.numeric(ncol(matrixRaw))
        rnum       <- as.numeric(nrow(matrixRaw))
        nnz        <- nnzero(matrixRaw)
        actionMessage(sprintf("%d x %d matrix, %.2f%% non-zero",
                              rnum, cnum, 100 * nnz / (cnum * rnum)),
                      prefix = "  ")
        rv
    },

    .readMatrixLOL = function () {
        ## Supporting an in-house flat file format for ragged list
        ## storage. This is really, really slow compared to reading
        ## from a .mtx file. For a 777 ID x 1335 list matrix, the
        ## times were:
        
        ## .mtx file : 910 msec (includes 240 msec of metadata parsing)
        ##  LOL file : 13.5 min, using SetFisherGrowableTriple (~ 900x slower)
        ##  LOL file : 20 sec, using identical SFGT code ripped out and
        ##             pasted inline with code (current state of function)
        
        ## The included Perl script listOfList_to_MatrixMarket.pl will
        ## convert the same file to MatrixMarket format in ~ 2 seconds
        
        dateMessage(paste("Reading List-of-Lists", colorize(file,"white")))
        fh     <- file(file, open = "r")
        block  <- 1000           # R Inferno trick - pre-expand vectors
        icnt   <- 0L             # Number of distinct rows
        inames <- integer(block) # Lookup for iName -> row number
        ivals  <- integer(block) # Triple component, i values
        jcnt   <- 0L             # Number of distinct columns
        jnames <- integer(block) # Lookup for iName -> col number
        jvals  <- integer(block) # Triple component, j values
        xvals  <- integer(block) # Triple component, x values
        xcnt   <- 0L             # Number of triples recorded

        jv <- 0 # Current list name index
        x  <- 0 # Current order in current list
        while (length(line <- readLines(fh, n = 1, warn = FALSE)) > 0) {
            if (line == "") next # Blank line
            if (grepl('^#', line)) {
                ## Comment line
                nameDat <- CatMisc::parenRegExp("^#\\s+LIST\\s+-\\s+(.+?)\\s*$", line)
                nm <- nameDat[1]
                if (!is.na(nm)) {
                    ## New list name
                    jv <- jnames[nm]
                    if (!is.na(jv)) {
                        ## Every LIST entry should be unique
                        err("List-of-List has multiple LIST entries with same name", nm, fatal = TRUE)
                    }

                    jcnt <- jcnt + 1L
                    ## Expand the vector if we have filled it up
                    if (jcnt > length(jnames))
                        jnames <- c(jnames, integer(block))
                    jnames[jcnt]        <- jcnt
                    names(jnames)[jcnt] <- nm
                    jv <- jcnt
                    x  <- 0
                }
                next
            }
            if (!grepl('^\\S', line)) next;
            ## Presume to be an id; Get the index value for this name
            iv <- inames[line]
            if (is.na(iv)) {
                ## We have not seen this entry yet, add it to the lookup
                icnt <- icnt + 1L
                ## Expand the vector if we have filled it up
                if (icnt > length(inames)) inames <- c(inames, integer(block))
                inames[icnt]        <- icnt
                names(inames)[icnt] <- line
                iv <- icnt
            }
            xcnt <- xcnt + 1L 
            if (xcnt > length(xvals)) {
                ## Need to expand the three triple vectors
                ivals <- c(ivals, integer(block))
                jvals <- c(jvals, integer(block))
                xvals <- c(xvals, integer(block))
            }
            ## Store the i,j,x triples
            ivals[ xcnt ] <- iv
            jvals[ xcnt ] <- jv
            xvals[ xcnt ] <- x  <- x + 1
        }
        close(fh)
        ## !!! Matrix wants zero-indexed columns!
        mat  <- new("dgTMatrix",
                    i   = ivals[1:xcnt] - 1L,
                    j   = jvals[1:xcnt] - 1L,
                    x   = as.numeric(xvals[1:xcnt]),
                    Dim = c(icnt,jcnt))
        rnDat         <- uniqueNames(names(inames[1:icnt]),
                                             "Row Names")
        cnDat         <- uniqueNames(names(jnames[1:jcnt]),
                                             "Col Names")
        rownames(mat) <- rnDat$names
        colnames(mat) <- cnDat$names

        list( matrix = mat,
             colChanges = cnDat$changes, rowChanges = cnDat$changes )
    },
    
    .readMatrixTXT = function ( listname = NA, header = FALSE ) {
        dateMessage(paste("Reading single list from text file", colorize(file,"white")))
        fh <- file(file, open = "r")
        if (header) {
            listname <- readLines(fh, n = 1, warn = FALSE)
        } else if (!CatMisc::is.something(listname)) {
            listname <- "MyList" # Eh, need something better
        }
        names <- character()
        while (length(line <- readLines(fh, n = 1, warn = FALSE)) > 0) {
            if (line == "" || grepl("^#", line)) next # Blank line or comment
            if (grepl('\\s', line)) {
                ## Normalize whitespace to single spaces
                line <- gsub('\\s+', ' ', line)
                ## Leading and trailing whitespace removal
                line <- gsub(' $', '', gsub('^ ', '', line))
            }
            names <- c(names, line)
        }
        close(fh)
        len  <- as.integer(length(names))
        len0 <- len - 1
        ## The constructor is VERY picky about what atomic classes are used here
        mat  <- new("dgTMatrix", i = 0:len0, j = rep(0L, len),
                    x = as.numeric(1:len), Dim = c(len, 1L))
        
        rnDat         <- uniqueNames(names, "Row Names")
        cnDat         <- uniqueNames(listname, "List Name")
        rownames(mat) <- rnDat$names
        colnames(mat) <- cnDat$names
        list( matrix = mat,
             colChanges = cnDat$changes, rowChanges = cnDat$changes )
    },
    
    .readMatrixGMT = function () {
        dateMessage(paste("Reading GMT set file", colorize(file,"white")))
        fh <- file(file, open = "r")
        names <- character()
        descr <- character()
        data  <- list()
        lnum  <- 0
        while (length(line <- readLines(fh, n = 1, warn = FALSE)) > 0) {
            row   <- unlist(strsplit(line, "\t"))
            names <- c(names, row[1])
            descr <- c(descr, row[2])
            lnum  <- lnum + 1
            data[[lnum]] <- row[ -(1:2) ]
        }
        close(fh)
        matrixFromLists(data, listNames=names,
                        meta=list(Description=setNames(descr,names)) )
    },

    .readMatrixMTX = function () {
        
        ## Parse row and column labels out of the comments in a
        ## MatrixMarket file. First read the file as a sparse matrix
        dateMessage(c("Reading Matrix file", colorize(file,"white")))
        mat <- readMM(file)
        rv  <- list()
        ## We will be taking a product, need to turn row/col counts to
        ## numeric to prevent "NAs produced by integer overflow"
        ## errors when dealing with large matrices
        nr  <- as.numeric(nrow(mat))
        nc  <- as.numeric(ncol(mat))
        msg <- sprintf("%d x %d, %.3g%% populated", nr, nc,
                       100 * nnzero(mat) / (nr * nc))
        ## Check to see if the file includes metadata
        ## First, find how many rows are comments
        lastCom <- max(grep("^%", readLines(file)))
        allCom  <- readLines(file, warn = FALSE, n = lastCom)
        ## It is possible there are comments later in the file? That
        ## is, we may have some non-comments in allCom
        comChk  <- grep("^%", allCom)
        if (length(comChk) < lastCom) allCom <- allCom[ comChk ]
        ## Remove leading comment token and leading whitespace
        allCom <- gsub('^%+\\s*', "", allCom)
        ## Is a separator defined?
        sep    <- grep("^Separator\\s+'.+'", allCom)
        if (length(sep) == 0) {
            sep <- NA
        } else if (length(sep) > 1) {
            err(c("Multiple separators defined in matrix file",
                  file, sep), fatal = TRUE)
        } else {
            sep <- CatMisc::parenRegExp("Separator\\s+\\'([^\\']+)\\'", allCom[sep])
        }
        fac    <- grep("^LEVELS\\s+\\[.+\\]", allCom)
        if (length(fac) == 0) {
            fac <- NA
        } else if (length(fac) > 1) {
            err(c("Multiple factor levels defined in matrix file",
                  file, fac), fatal = TRUE)
        } else {
            ## Not sure why "[^]]" works - I should have to escape the
            ## internal bracket (eg "[^\\]]" but that does NOT work.
            facDat <- CatMisc::parenRegExp("LEVELS\\s+\\[([^]]+)\\]\\[([^]]+)\\]",
                               allCom[fac])
            rv$levels <- unlist(base::strsplit(facDat[2], facDat[1]))
        }
        
        ## Find the row and column header lines
        rcPos <- grep("^(Row|Col) Name", allCom)
        rcLen <- length(rcPos)
        ## Get any defaults that have been set
        def <- grep("^DEFAULT\\s+", allCom)
        params <- list()
        for (defl in def) {
            kv <- CatMisc::parenRegExp("DEFAULT\\s+(\\S+)\\s+(\\S.*?)\\s*$", allCom[defl])
            params[[ tolower(kv[1]) ]] <- kv[2]
        }
        rv$params <- params
        parseRow <- function (x, sep = NA) {
            ## strsplit() on spaces discards trailing spaces, which
            ## causes problems for a space-containing metadata token
            ## (eg the default ' :: ')when the last metadata column is
            ## empty. So use CatMisc::parenRegExp to get the index number instead:

            indDat <- CatMisc::parenRegExp("^(\\d+)\\s+(.+)", x)
            if (is.na(sep)) {
                ## No metadata separator defined
                indDat
            } else {
                ## Split the data section on the separator
                c(indDat[1], unlist(base::strsplit(indDat[2], sep)) )
            }

            
            
            ## Break on spaces to get index num at front
            ## This is NOT vectorized, due to unlist()
            ## s <- unlist(base::strsplit(x, "\\s+"))
            ## Reconstitute the other data. Note if the metadata had
            ## space runs in it, these will be collapsed to single
            ## spaces.
            ## m <- paste(s[-1], collapse = " ")
            ## if (!is.na(sep)) m <- unlist(base::strsplit(m, sep))
            ## c(s[1], m)
        }

        padRow <- function (x, len = 0) {
            ## Make sure a vector is at least len long, padding with NAs
            if (length(x) < len) x[len] <- NA
            x
        }
        rChng <- NULL
        cChng <- NULL
        metadata <- data.table::data.table( id = character(), key = "id" )
        dimNames <- list( Row = NULL, Col = NULL )
        dimChngs <- list( Row = NULL, Col = NULL )
        for (i in seq_len(rcLen)) {
            ## Is this specifying Row or Col, and what are the headers?
            targDat <- CatMisc::parenRegExp("^(\\S+)\\s+(.+?)\\s*$", allCom[rcPos[i]])
            what    <- targDat[1] # Row / Col
            dateMessage(sprintf("Parsing %s names", what), prefix="  ")
            meta <- targDat[2]
            if (!is.na(sep))
                meta <- base::strsplit(meta, sep, fixed = TRUE)[[1]]
            meta[1] <- "id" # Generally will be "Name", but normalize to "id"
            ## What are the first and last coordinate rows?
            s <- rcPos[i] + 1
            e <- ifelse(i == rcLen, length(allCom), rcPos[i+1] - 1)
            rows <- allCom[s:e]
            ## Remove non-coordinate rows
            isCoord <- grep("^\\d+\\s+\\S", rows)
            ## Make a (possibly) ragged list
            ragged  <- lapply(rows[ isCoord ], parseRow, sep = sep)
            ## Find the maximum number of columns
            maxCol  <- max(unlist(lapply(ragged, length)))
            ## Pad and convert to character matrix
            cmat    <- base::matrix(unlist(lapply(ragged, padRow, maxCol)),
                                    ncol = maxCol, byrow = TRUE)
            tmpIndName     <- ".index."
            colnames(cmat) <- uniqueNames(
                c(tmpIndName, meta), paste("Metadata headers for",what))$names
            rownames(cmat) <- uniqueNames(cmat[, "id" ])$names
            ## The first column are indices in the sparse matrix
            inds    <- as.integer(cmat[, 1])
            ## The second column holds the names
            names <- character()
            names[ inds ] <- cmat[, 2]
            mnn <- uniqueNames(names,paste(what, "Names"))
            dimNames[[ what ]] <- mnn$names
            dimChngs[[ what ]] <- mnn$changes
            mlen <- length(meta)
            if (mlen > 1) {
                ## Reliably building the DT proved
                ## non-intuitive. Using the matrix directly with
                ## as.data.table() seems to be the most reliable
                ## mechanism.

                ## https://stackoverflow.com/a/10235618
                ## Convert the matrix to a DT, leaving out the index column
                tmp <- data.table::as.data.table( cmat[, -1, drop = FALSE],
                                                 key = "id")
                rownames(tmp) <- uniqueNames(cmat[, "id" ])$names
                ## I never did find a way to add *ROWS* by reference in DTs...
                metadata <- data.table::rbindlist(list(metadata, tmp),
                                                  fill = TRUE)
                data.table::setkey(metadata, "id") ## Merged DTs don't carry over key
            }
        }
        ## Set the dimension names
        rdn <- params$rowdim
        names(dimNames)[1] <- ifelse(CatMisc::is.something(rdn), rdn, "")
        cdn <- params$coldim
        names(dimNames)[2] <- ifelse(CatMisc::is.something(cdn), cdn, "")
        dimnames(mat) <- dimNames

        ## Now check if there are also sidecar metadata files. We will
        ## expect them to have the same name as the MatrixMarket file,
        ## but be additionally suffixed with
        ## "-metadata<any-other-characters>"
        fileDir   <- dirname(file)
        baseFile  <- basename(file)
        metafiles <- list.files(fileDir, ignore.case=T,
                                pattern = sprintf("^%s-metadata", baseFile))
        for (mf in metafiles) {
            ## Explicit "sidecar" metadata are present
            ## Read in as temporary data.table
            dateMessage(c("Reading metadata file -", colorize(mf,"white")),
                        prefix="  " )
            scDT <- data.table::fread(mf)
            ## Find the "id" column, or the first one
            scNm <- names(scDT)
            scIdCol <- which(scNm == "id")
            if (length(scIdCol) == 0) {
                ## No column named 'id'. *Presume* it is the first
                ## one and rename it.
                setnames(scDT, scNm[1], "id")
            }
            data.table::setkey(scDT, "id")
            ## Set the key to "id", and then merge into metadata
            metadata <- data.table::rbindlist(list(metadata, scDT), fill = TRUE,
                                  use.names = TRUE)
            data.table::setkey(metadata, "id")
        }
        rv$matrix     <- mat
        rv$metadata   <- metadata
        rv$colChanges <- dimChngs$Col
        rv$rowChanges <- dimChngs$Row
        rv
     },
    
     metadata = function ( id = NULL, key = NULL) {
        if (!CatMisc::is.something(key)) {
            ## No column specified
            if (!CatMisc::is.something(id)) return( NULL )
            ## Return a subset of the data table:
            rv <- matrixMD[ id, ]
            data.table::setkey(rv, "id")
            return( rv )
        } else if (!CatMisc::is.something(id)) {
            ## Return a named vector for the metadata column
            allIds <- as.character(matrixMD[,("id"),with=FALSE][[1]])
            return(setNames( matrixMD[allIds,key,with=FALSE][[1]], allIds))
        }
        ## Else we are asking for specific metadata for specific IDs
        matrixMD[ id, key, with = FALSE ]
    },

    getRow = function ( x = NA, format = "vector", sort = NA,
        min = NA, max = NA) {
        mat <- matrix()
        ## Get the rows we are interested in:
        rv  <- mat[ x, , drop = FALSE]
        ## Apply (by "masking" to zero) min/max limits, if requested
        if (!is.na(min)) rv[ rv < min ] <- 0
        if (!is.na(max)) rv[ rv > max ] <- 0
        ## Eliminate all zero cols: https://stackoverflow.com/a/6632287
        rv  <- rv[ , colSums(abs(rv)) != 0, drop = FALSE ]
        
        
        fmt <- tolower(substr(format, 1, 1))
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
            ## Just return the matrix
            rv
        }
    },

    getCol = function ( x = NA, format = "vector", sort = NA,
        min = NA, max = NA) {
        mat <- matrix()
        ## Get the columns we are interested in:
        rv  <- mat[ , x, drop = FALSE ]
        ## Apply (by "masking" to zero) min/max limits, if requested
        if (!is.na(min)) rv[ rv < min ] <- 0
        if (!is.na(max)) rv[ rv > max ] <- 0
        ## Eliminate rows of all zeros: https://stackoverflow.com/a/6632287
        rv  <- rv[ rowSums(abs(rv)) != 0, , drop = FALSE ]
        
        fmt <- tolower(substr(format, 1, 1))
        if (fmt == "v") {
            ## vector
            ## Just return a vector of matching column names
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

    matrixText = function ( pad = "", useObj = NULL, fallbackVar = NULL,
                          compact = FALSE, color=NULL ) {
        if (is.null(color)) color <- useColor() # Use EventLogger setting
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        ## Variable name for use in sample methods
        objName <- .self$.selfVarName("myMatrix", fallbackVar)
        msg <- doCol(sprintf("%s sparse matrix\n", class(.self)[1]), "blue")
        ## Even if we read the file from RDS, we will base the
        ## reported age on the original data file (if available)
        src <- ifelse(file.exists(file), file, sprintf("%s.rds", file))
        age <- sprintf("%.1f",difftime(Sys.time(),file.mtime(src),units="days"))
        msg <- sprintf("%s  Source: %s%s [%s days old]\n",
                       msg, doCol(file, "white"),
                       ifelse(fromRDS, doCol('.rds', "purple"),""),
                       doCol(age,"red"))
        if (is.factor() && !compact) {
            ## Report the factor levels
            indent <- "\n      ";
            lvl <- paste(doCol(strwrap(paste(lvlVal, collapse = ', '),
                                          width = 0.7 * getOption("width")),
                                  "purple"), collapse = indent)
            msg <- sprintf("%s    %s%s%s\n", msg, doCol(
                "Values are factors:", "yellow"), indent, lvl)
        }
        mNm <- attr(matrixRaw, "matrixName")
        name <- param("name")
        if (CatMisc::is.something(name))
            msg <- sprintf("%s    Name: \"%s\"\n", msg, doCol(name, "white"))
        desc <- param("description")
        if (CatMisc::is.something(desc) && !compact)
            msg <- sprintf("%s    Desc: %s\n", msg, doCol(desc, "white"))
        
        ## If an external "used" matrix is not provided, use internal field:
        if (!CatMisc::is.def(useObj)) useObj <- matrixUse
        dimNames <- names(dimnames(matrixRaw))
        mat <- if (is.null(useObj)) { matrixRaw } else { useObj }
        nr  <- nrow(mat)
        nc  <- ncol(mat)
        nz  <- nnzero(mat)
        msg <- sprintf("%s  %8d %s", msg, nr, ifelse(CatMisc::is.something(
            dimNames[1]),dimNames[1],"rows"))
        rN  <- rownames(mat)
        if (CatMisc::is.def(rN))
            msg <- sprintf("%s eg: %s", msg, doCol(substr(paste(rN[1:pmin(3,length(rN))], collapse = ', '), 1, 60), "cyan"))
        
        msg <- sprintf("%s\n  %8d %s", msg, nc, ifelse(CatMisc::is.something(
            dimNames[2]),dimNames[2],"cols"))
        cN  <- colnames(mat)
        if (CatMisc::is.def(cN))
        msg <- sprintf("%s eg: %s", msg, doCol(substr(paste(cN[1:pmin(3,length(cN))], collapse = ', '), 1, 60), "cyan"))

        perPop <- sprintf("%.3g%%", 100 * nz / (nr * nc))
        msg <- sprintf("%s\n  %d non-zero cells (%s)", msg, nz,
                       doCol(perPop, "red"))
        msg <- sprintf("%s\n", msg)
        ## Note if we had to remap one or more R/C names
        changed <- character()
        if (CatMisc::is.def(rowChanges)) changed <- c(changed, doCol(
            paste(length(rowChanges), "Rows ($rowChanges)"), "yellow"))
        if (CatMisc::is.def(colChanges)) changed <- c(changed, doCol(
            paste(length(colChanges), "Cols ($colChanges)"), "yellow"))

        if (length(changed) != 0) msg <- sprintf(
                      "%s    %s names needed to be altered.\n",
                      msg, paste(changed, collapse = " and "))
        
        if (!is.null(useObj)) {
            ## Show stats for raw matrix
            ## Numeric conversion to prevent integer overflow on product
            rnr <- as.numeric(nrow(matrixRaw))
            rnc <- as.numeric(ncol(matrixRaw))
            rnz <- nnzero(matrixRaw)
            rpp <- sprintf("%.3g%%", 100 * rnz / (rnr * rnc))
            ## Note the percent change for Use/Raw
            rcD <- doCol(sprintf("%.1f%%x%.1f%%", 100 * nr/rnr,
                                    100 * nc / rnc), "yellow")
            rzD <- doCol(sprintf("%.1f%%", 100 * nz / rnz), "yellow")
            msg <- sprintf("%s  %s : %d x %d (%s), %d non-zero (%s; %s survive)\n",
                           msg, doCol("[Raw Matrix]", "red"), rnr, rnc, rcD,
                           rnz, doCol(rpp, "red"), rzD)
        }
        mdRow <- ifelse(is.null(matrixMD), 0, nrow(matrixMD))
        if (mdRow != 0 && !compact) {
            mcols <- names(matrixMD)[-1]
            msg   <- paste(c(msg, doCol(sprintf("  %d IDs have metadata assigned in up to %d keys, eg:\n", mdRow, length(mcols)), "blue"), collapse = ""))
            for (mcol in mcols) {
                ## Get a non-NA sample ID with this annotation
                notNA <- !is.na(matrixMD[[ mcol]])
                ## OMG ITS SO PAINFUL JUST TO GET A VECTOR
                exid  <- as.character(matrixMD[ notNA,1,with = FALSE ][[1]][1])
                exval <- matrixMD[ exid, mcol, with = FALSE][[1]]
                mline <- sprintf('    %s$metadata(key="%s", id="%s") > "%s"\n',
                                 doCol(objName, "white"),
                                 doCol(mcol,"red"),
                                 doCol(exid, "cyan"),
                                 doCol(exval, "magenta"))
                msg <- paste(c(msg, mline), collapse ="")
            }
        }
        ## Left pad if requested
        if (CatMisc::is.def(pad)) msg <-
            paste(c(sprintf("%s%s", pad, unlist(base::strsplit(msg, "\n"))),
                    "\n"), collapse = "\n")
        ## Ending up with extraneous newlines at end
        sub("\n+$", "\n", msg)
    }
    
)


#' Files to MatrixMarket
#'
#' Converts a collection of files to a single MatrixMarket file
#'
#' @param files Required, a vector of file paths
#' @param output Default "SparseMatrix.mtx", the filename of the
#'     MatrixMarket file being generated
#' @param filter Optional function used to filter set elements. It
#'     should take as input a character vector and return a character
#'     vector representing the desired (kept) elements.
#' @param rename Optional function used to generate list names. It
#'     should take as input the file path (single string) and return a
#'     single string that will be used as a name for that list
#'
#' @examples
#' 
#' \dontrun{
#'   # Find all files names "geneList#.txt" and extact Ensembl Gene IDs
#'   # from them into a single MTX file
#'   filesToMTX(files=list.files(pattern='geneList[0-9]+.txt'),
#'              filter=function(x) x[ grepl("^ENSG\\d+$", x) ],
#'              output="myEnsemblLists.mtx")
#' }
#' 
#' @export

filesToMTX <- function(files, output="SparseMatrix.mtx",
                       filter=NULL, rename=NULL) {
    ## Convert a set of text files to a MatrixMarket file
    cn <- character()
    data <- list()
    for (file in files) {
        n <- if (is.null(rename)) { file } else { rename(file) }
        if (n %in% cn) {
            message("File ", n, " requested multiple times")
            next
        }
        cn <- c(cn, n)
        l  <- utils::read.table(file, stringsAsFactors=F)
        dl <- length(data)
        data[[ dl+1 ]] <- if (is.null(filter)) { l[[1]] } else {filter(l[[1]])}
    }
    allEntries <- unlist(data)
    allIds     <- unique(allEntries)
    tot        <- length(allEntries)
    nr         <- length(allIds)
    nc         <- length(cn)
    
    fh <- file(output, "w")
    writeLines(c("%%MatrixMarket matrix coordinate real general",
                 "% Row Name",
                 sprintf("%% %d %s",seq_len(nr),allIds),
                 "% Col Name",
                 sprintf("%% %d %s",seq_len(nc),cn),
                 "% Matrix triples : Row Col Score",
                 sprintf("  %d %d %d", nr, nc, tot)),fh)
    for (col in seq_len(nc)) {
        rows <- match(data[[col]], allIds)
        numNA <- sum(is.na(rows))
        if (numNA > 0) {
            message(c(cn[col], "Some entries could not be match to full set",
                     "MatrixMarket file will likely fail to load"))
            rows <- rows[ !is.na(rows) ]
        }
        ## The 'score' will be the order the ID appeared in the original list
        score <- match(allIds[rows], data[[col]])
        if (length(score) != length(rows)) {
            message("Somehow managed to mess up query rank calculation")
        }
        writeLines(sprintf("%d %d %d", rows, col, score), fh)
    }
    message("MatrixMarket file written to ", output)
    close(fh)
    invisible(output)
}

#' Take Lowest Thing
#'
#' From a list of strings, take the lowest/smallest thing
#'
#' @details
#'
#' This is a desperation function and should only be used when you
#' have a character vector and absolutely must reduce it to
#' "one thing" but you have no real mechanism to do so. It was
#' designed to address requests to reduce a list of gene accessions
#' (eg LOC13992, LOC93) to "one gene", but without any other guiding
#' context. The function will extract the left-most uninterupted
#' integer that it can find, and pick the smallest. The rationale is
#' that "LOC93" was probably annotated earlier (longer ago) than
#' "LOC13200423", and as such is probably the "more common" /
#' "better annotated" / "more popular" of the two.
#'
#' @param x A character vector of strings to pick from
#'
#' @return A single string
#'
#' @examples
#'
#' x <- c("LOC828221", "LOC1234", "LOC39", "HUH?", "LOC39-BETA", "LOC99.243-X")
#' takeLowestThing(x)
#' 
#' @export

takeLowestThing <- function (x) {
    ## Grab the 'first' full integer out of each string:
    nums <- as.integer(gsub('[^0-9].*', '', gsub('^[^0-9]+', '', x)))
    ## Which one(s) are the smallest?
    inds <- which(nums == min(nums, na.rm=TRUE))
    ## Sort (alphabetically) and return the 'smallest'
    sort(x[inds])[1]
}

#' Unique Names
#'
#' Wrapper for make.names and make.unique, with reporting of changes
#'
#' @param names Required, a character vector of the names to process
#' @param rpt Default NULL. If not NULL, will report to STDERR any
#'     changes that were required. The value of rpt will be used as a
#'     noun in the message (ie if rpt='gene' the message will be of
#'     the form "6 genes required alteration: "...)
#' @param valid Default FALSE, whch will utilize make.unique. If TRUE,
#'     then make.names( unique=TRUE) will be used instead
#'
#' @return
#'
#' A list with two components: "names", the vector of (possibly)
#' transformed names, and "changes", a named vector where the names
#' are the new names and the values are the original ones (will only
#' include altered names)
#'
#' @examples
#'
#' foo <- c("Apple","Banana","Cherry","Apple")
#' uniqueNames(foo, rpt="IDs")
#'
#' @importFrom crayon bgCyan
#' 
#' @export

uniqueNames <- function(names = character(), rpt = NULL,
                         valid = FALSE ) {
    ## Normalize names for rows/columns, and also record (and
    ## optionally report) any differences
    goodNames <- NULL
    if (valid) {
        ## Also make the names valid
        goodNames <- make.names( names, unique = TRUE )
    } else {
        goodNames <- make.unique( names )
    }
    changes   <- NULL
    if (!identical(names, goodNames)) {
        ## Some changes were made. Note them as a (good)named subset
        diff    <- goodNames != names
        changes <- setNames(names[diff], goodNames[diff])
        if (!is.null(rpt)) {
            ## Note the changes to STDERR
            num <- sum(diff)
            msg <- crayon::bgCyan(paste(num,rpt,"required alteration"))
            if (num <= 20) msg <- c(msg, vapply(seq_len(length(changes)),
                   function (i) {
                       sprintf("  '%s' -> '%s'",changes[i],
                               names(changes)[i])
                   }, ""))
            message(msg, collapse = "\n")
        }
    }
    list( names = goodNames, changes = changes )
}

#' Matrix from Lists
#'
#' Convert a list of vectors into a sparse matrix consumable by AnnotatedMatrix
#'
#' @details
#'
#' Sparse matrices (eg dgTMatrix) are constructed from lists of
#' triples; row number (i), column number (j) and value of that cell
#' (x). This function generates such triples by providing one or more
#' lists of "members" (strings)
#'
#' @param data Required, a list of vectors. Each list element is
#'     considered to be a row, and the members will become (in
#'     aggreagate) the columns
#' @param meta Default NULL. Optional list of named character vectors
#'     representing metadata assignements. The list names represent
#'     the metadata names (eg "Description", "Notes", whatever). The
#'     vector names represent the thing being annotated (a row or
#'     column name), and the values hold the actual
#'     assignment. Metadata is "mixed" in that there are not separate
#'     tables for columns and rows.
#' @param listNames Default NULL, in which case the listNames
#'     (rownames) will be taken as names(data). The option to provide
#'     the names independently is included in case some names are not
#'     unique. They will still be made so in the final matrix
#' @param val Default 1. The value to be assigned to non-zero cells in
#'     the matrix.
#'
#' @return
#'
#' A list with the following members:
#'
#' \itemize{
#'
#'   \item mat - The sparse matrix
#'
#'   \item metadata - data.table of any metadata assignments
#'
#'   \item colChanges - a named character vector. The names are the
#'   column names as used in the matrix, the values are the names as
#'   originally provided. Will only contain names that had to be
#'   altered to avoid uniqueness conflicts.
#'
#'   \item rowChanges - as above, but for rows
#'
#' }
#'
#' This structure can be read by internal methods used by
#' AnnotatedMatrix to parse flat files
#'
#' @examples
#'
#' ingredients <- list(Potato=c("Mashed","Fried","Baked"),
#'                     Onion=c("Fried"),
#'                     Apple=c("Raw","Baked"))
#' 
#' meta <- list(Type=setNames(c("Vegetable","Fruit","Vegetable"),
#'                            c("Potato","Apple","Onion")),
#'              Tool=setNames(c("Tray","Pan","Pot",NA),
#'                            c("Baked","Fried","Mashed","Raw")),
#'              Appliance=setNames(c("Oven","Stove","Stove",NA),
#'                                 c("Baked","Fried","Mashed","Raw")))
#'                              
#' matrixFromLists(ingredients, meta)
#'
#' @importClassesFrom Matrix dgTMatrix
#' @importFrom data.table as.data.table setkeyv
#' @importFrom stats setNames
#' 
#' @export

matrixFromLists <- function(data, meta=NULL, listNames=NULL, val=1) {
    ## Take list (row) names from data if not explicitly provided:
    if (is.null(listNames)) listNames <- names(data)
    rnDat  <- uniqueNames(listNames, "Row Names")
    icnt   <- length(rnDat$names)
    ## All unique column names:
    memNames <- unique(unlist(data))
    cnDat    <- uniqueNames(memNames, "Col Names")
    jcnt     <- length(cnDat$names)
    ## How many members are in each list?
    counts   <- unlist(lapply(data, length))
    ## Total number of assingments:
    xcnt     <- sum(counts)
    ## Pre-allocate the i/j index vectors:
    ivals    <- integer(xcnt)
    jvals    <- integer(xcnt)
    done     <- 0
    for (l in seq_len(icnt)) {
        ## Get the set of members from List #l :
        mems <- data[[ l ]]
        mlen <- length(mems)
        inds <- seq(done+1, done+mlen) # Indices for i/j vectors
        ivals[ inds ] <- rep(l, mlen)  # Set i coordinate
        jvals[ inds ] <- match(mems, memNames) # Set j coordinate
        done <- done + mlen
    }
    ## !!! Matrix wants zero-indexed columns!
    mat  <- new("dgTMatrix",
                i   = ivals[1:xcnt] - 1L,
                j   = jvals[1:xcnt] - 1L,
                x   = rep(val, xcnt),
                Dim = c(icnt,jcnt))
    rownames(mat) <- rnDat$names
    colnames(mat) <- cnDat$names

    metadata <- NULL
    if (!is.null(meta)) {
        ## Get all the unique IDs for which there is metadata
        ids   <- unique(unlist(lapply(meta, names)))
        ## Make columns
        mList    <- lapply(meta, function(x) x[ids])
        mList$id <- ids
        metadata <- data.table::as.data.table( mList,  key = "id")
        data.table::setkeyv(metadata, "id")
    }

    list(matrix = mat, metadata=metadata,
         colChanges = cnDat$changes, rowChanges = cnDat$changes )
}

