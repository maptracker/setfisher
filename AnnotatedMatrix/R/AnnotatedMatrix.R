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
#' @importFrom data.table data.table as.data.table
#'      setkeyv
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
    
    matObj = function( raw=FALSE, transpose=FALSE) {
        "\\preformatted{
Retrieves the underlying Matrix for this object. Parameters:
      raw - Default FALSE, in which case the filtered Matrix (held in field
            'matrixUse') will be returned, if it is available. If not
            available, or if raw is TRUE, then the raw (as loaded from file)
            Matrix will be returned.
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

    rNames = function( raw=FALSE ) {
        "\\preformatted{
Returns the row names of the Matrix
      raw - Default FALSE, in which case the filtered Matrix (held in field
            'matrixUse') will be used, if available - otherwise the original
            raw matrix will be used.
}"
        rownames(matObj(raw))
    },

    cNames = function( raw=FALSE ) {
        "\\preformatted{
Returns the column names of the Matrix
      raw - Default FALSE, in which case the filtered Matrix (held in field
            'matrixUse') will be used, if available - otherwise the original
            raw matrix will be used.
}"
        colnames(matObj(raw))
    },

    reset = function( asFactor=FALSE ) {
        "\\preformatted{
Reset any filters that were applied - the 'used' matrix will be the
original 'raw' one
}"
        matrixUse <<- NULL
        filterLog <<- data.table(id = character(), key = "id")
        invisible(NA)
    },

    filterByScore = function( min=NA, max=NA, filterEmpty=FALSE, reason=NA ) {
        "\\preformatted{
Apply filters to the current matrix to zero-out cells failing thresholds.
Returns the number of cells zeroed (filtered) out.
        min - Minimum allowed value. Cells below this will be set to zero
        max - Maximum allowed value. Cells above it will be set to zero
 filterEmpty - Default FALSE; If true, then the matrix will be 'shrunk' to
              remove rows and columns that are only zeros
     reason - Default NA; If specified, a text value that will be added to
              the $filterLog under the 'reason' column
}"
        obj <- matObj()
        rv  <- 0
        if (filterEmpty) removeEmpty("Empty rows and cols before score filter")
        if (is.something(min)) {
            ## Zero out entries that fall below min
            fail <- obj@x < min
            numZ <- sum(fail)
            if (numZ > 0) {
                ## At least some cells were zeroed out
                obj@x[ fail ] <- 0
                ## obj <- drop0(obj)
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
                    obj <- matObj()
                }
            }
        }
        if (is.something(max)) {
            ## Zero out entries that are above max
            fail <- obj@x > max
            numZ <- sum(fail)
            if (numZ > 0) {
                ## At least some cells were zeroed out
                obj@x[ fail ] <- 0
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
                    obj <- matObj()
                }
            }
        }
        if (rv != 0 && !filterEmpty) matrixUse <<- obj
        invisible(rv)
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
        obj <- matObj()
        rv  <- 0
        xVal <- integer()
        lvls <- levels()
        nlvl <- length(lvls)
        if (is.numeric(x)) {
            xVal <- if (!is.integer( x )) {
                xVal <- as.integer( x )
                if (!all.equal(x, xVal)) {
                    err(c("Integer conversion in filterByFactorLevel() appears to have altered your request"))
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


        ## Working with the @x value vector for the sparse
        ## matrix. Determine which values are 'unwanted'
        fail <- is.element(obj@x, xVal)
        if (keep) fail <- !fail
        rv <- sum( fail )
        if (rv) {
            ## Zero-out the failing cells
            if (filterEmpty) removeEmpty("Empty rows and cols prior to factor filter")
            obj@x[ fail ] <- 0
            ## Even if the user asked to keep a set of levels, the set
            ## being excluded might be smaller (and vice
            ## versa). Determine what the most compact representation
            ## of the filter is
            k  <- keep
            lnames <- if (length(xVal) * 2 > nlvl) {
                ## Briefer to refer to the un-specified set
                k <- !k
                lvls[ -xVal ]
            } else {
                lvls[ xVal ]
            }
            ## The filter text describes failing cells
            testTxt <- sprintf("x %sassigned to %s", if (k) {"not "} else {""},
                               paste(lnames, collapse=', '))
            numTxt  <- paste(rv, "Cells")
            filterDetails(id=testTxt, type="Factor", metric=numTxt,
                          reason=reason)
            if (filterEmpty) {
                ## Strip empty rows and columns
                matrixUse <<- obj
                if (!is.na(reason[1])) testTxt <- paste(testTxt, reason[1])
                removeEmpty(testTxt)
                obj <- matObj()
            }
        }
        if (rv != 0 && !filterEmpty) matrixUse <<- obj
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
        if (length(toss) != 0) {
            ## Some columns have been removed
            filterDetails(id=toss, type="Row", metric="AllZero", reason=reason)
            matrixUse <<- obj[ !isEmpty, ]
        }
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

    populatedColumns = function(obj=NULL, ...) {
        "\\preformatted{
Return a logical vector indicating which columns have at least one non-zero
 cell
        obj - Default NULL, which will recover the matrix from matObj(),
              passing ... as well (so you can select the raw matrix if
              desired)
}"
        cCounts( obj=obj, ...) != 0
    },

    removeEmptyColumns = function(reason=NA) {
        "\\preformatted{
Remove all empty columns (those that only contain zeros). Invisibly
returns a vector of removed IDs.
     reason - Default NA; If specified, a text value that will be added
              to the $filterLog under the 'reason' column
}"
        obj     <- matObj()
        isEmpty <- !populatedColumns()
        ## populatedColumns() returns a named logical vector. Take
        ## names from there
        toss    <- names(isEmpty)[ isEmpty ]
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

    map = function(input, via=NULL, format="data.frame",
                   ignore.case=TRUE, keep.best=FALSE,
                   column.func=max,
                   collapse=NULL, collapse.name=NULL, collapse.token=',',
                   collapse.score=NULL,
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
}"
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
        valueCollapseFunction <-
        if (!is.null(collapse.score)) {
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

        if (is.something(add.metadata)) {


### TODO: Need to manage cases where colIn is true, and Output is a
### token-separated string of IDs. That is, need to grab metadata in
### some fashion for these synthetic values. Should probably hold a
### temp list for each row of original values
            
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

        if (grepl('vec', format[1], ignore.case=TRUE)) {
            ## Vector format
            rv <- setNames(rv$Output, rv$Input)
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
        row <- data.table::data.table(id=as.character(id), metric=metric,
                                      type=type, reason=reason, key="id")
        ## Do not add any entries that are already recorded as
        ## filtered (only note the first exclusion)
        newR <- base::setdiff(row[["id"]], filterLog[["id"]])
        filterLog <<- extendDataTable(filterLog, row[newR])
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
            fname <- file
            is.gz <- CatMisc::parenRegExp('(.+)\\.gz$', fname)
            if (!is.na(is.gz[1])) fname <- is.gz[1]
            
            if (grepl('(mtx|matrixmarket)', format, ignore.case = TRUE) ||
                grepl('\\.mtx', fname, ignore.case = TRUE)) {
                rv <- .readMatrixMTX( ... )
            } else if (grepl('(txt|text)', format, ignore.case = TRUE) ||
                       any(dir.exists(file)) ||
                       grepl('\\.(txt|text|list)$',fname, ignore.case = TRUE)) {
                dateMessage(paste("Reading simple list from text file",
                                  colorize(file,"white")))
                rv <- parse_Text_file( file=file, ... )
            } else if (grepl('(lol)', format, ignore.case = TRUE) ||
                       grepl('\\.(inp)$', fname, ignore.case = TRUE)) {
                dateMessage(paste("Reading List-of-Lists",
                                  colorize(file,"white")))
                rv <- parse_ListOfLists_file( file )
            } else if (grepl('(gmt)', format, ignore.case = TRUE) ||
                       grepl('\\.(gmt)$', fname, ignore.case = TRUE)) {
                dateMessage(paste("Reading lists from GMT file",
                                  colorize(file,"white")))
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
        actionMessage(sprintf("%d x %d matrix, %.2f%% non-zero",
                              rnum, cnum, 100 * nnz / (cnum * rnum)),
                      prefix = "  ")
        rv
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
            cmat    <- matrix(unlist(lapply(ragged, padRow, maxCol)),
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
                metadata <- extendDataTable(metadata, tmp)
                ## Merged DTs don't carry over key:
                data.table::setkeyv(metadata, "id") 
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

        metadata      <- parseMetadataSidecar( file, metadata )
        rv$matrix     <- mat
        rv$metadata   <- metadata
        rv$colChanges <- dimChngs$Col
        rv$rowChanges <- dimChngs$Row
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
            unknown <- setdiff(key, metadata_keys())
            if (length(unknown) > 0) {
                key <- intersect(key, metadata_keys())
                if (verbose) err(c("Request for unknown metadata key(s):",
                                   unknown))
                if (length(key) == 0) {
                    key <- NULL
                    if (is.null(id)) return( rv )
                }
            }
            ## Make sure the id column is included
            useKey <- c("id", key)
        }
        ## data.table does *not* like numeric keys...
        if (!is.null(id)) id <- as.character(id)

        if (is.null(id)) {
            ## Key request only
            rv <- matrixMD[ , useKey, with = FALSE ]
            if (na.rm) {
                ## Remove null rows. Find NA values in columns:
                naCol <- sapply(key, function(x) is.na( rv[[ x ]] ) )
                ## And now find rows that are all() NAs:
                naRow <- apply(naCol, 1, all)
                ## ... and remove them:
                rv <- rv[ naRow, ]
            }
        } else if (is.null(key)) {
            ## ID request only
            rv <- matrixMD[ id, , with = FALSE ]
            if (na.rm) {
                ## Remove null columns. Find NA values in columns:
                cols   <- colnames(rv)
                naCol1 <- sapply(cols, function(x) is.na( rv[[ x ]] ) )
                ## And then columns that are all null:
                naCol2 <- apply(naCol1, 2, all)
                ## ... and remove them:
                rv <- rv[ , naCol2 ]
            }
        } else {
            ## Both
            rv <- matrixMD[ id, useKey, with = FALSE ]
            ## We will not remove either rows or columns, since both
            ## were explicitly requested
        }

        cns <- colnames(rv)
        if (drop && length(cns) <= 2) {
            names <- rv[[ "id" ]]
            if (length(cns) == 1) {
                ## Big batch of nothing
                rv <- setNames(rep(NA, length(names)), names)
            } else {
                ## What metadata column do we have?
                mCol <- setdiff(cns, "id")
                rv   <- setNames(rv[[ mCol[1] ]], names)
            }
        }
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
        rv  <- obj[ x, , drop = FALSE]
        ## Apply (by "masking" to zero) min/max limits, if requested
        if (!is.na(min)) rv[ rv < min ] <- 0
        if (!is.na(max)) rv[ rv > max ] <- 0
        ## Eliminate all zero cols: https://stackoverflow.com/a/6632287
        rv  <- rv[ , Matrix::colSums(abs(rv)) != 0, drop = FALSE ]
        
        
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
        rv  <- obj[ , x, drop = FALSE ]
        ## Apply (by "masking" to zero) min/max limits, if requested
        if (!is.na(min)) rv[ rv < min ] <- 0
        if (!is.na(max)) rv[ rv > max ] <- 0
        ## Eliminate rows of all zeros: https://stackoverflow.com/a/6632287
        rv  <- rv[ Matrix::rowSums(abs(rv)) != 0, , drop = FALSE ]
        
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
        nr  <- as.numeric(nrow(mat))
        nc  <- as.numeric(ncol(mat))
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
            mcols <- metadata_keys()
            msg   <- paste(c(msg, doCol(sprintf("  %d IDs have metadata assigned in up to %d keys, eg:\n", mdRow, length(mcols)), "blue"), collapse = ""))
            for (mcol in mcols) {
                ## Get a non-NA sample ID with this annotation
                notNA <- !is.na(matrixMD[[ mcol]])
                rowNum <- which(notNA)[1]
                ## OMG ITS SO PAINFUL JUST TO GET A VECTOR
                exid  <- as.character(matrixMD[ rowNum, "id", with=FALSE][[1]])
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

### - Functions - ###

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
#' @param verbose Default TRUE, which will emit a warning if any names
#'     needed to be changed
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
                         valid = FALSE, verbose=TRUE ) {
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
        if (!is.null(rpt) && verbose) {
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

#' Make Temp File
#'
#' Internal method to help manage example files in extdata/
#'
#' @details
#'
#' When matrix files are loaded, a .rds version is generated to aid in
#' rapid loading in the future. These files are made alongside the
#' original flat file. This could cause problems if the file is from
#' inst/extdata, so this method takes a requested file and copies it
#' to a temporary location before loading.
#'
#' @param name Required, the basename of the file
#'
#' @examples
#'
#' tmpFile <- AnnotatedMatrix::.makeTempFile("Symbol-To-Gene.mtx")
#'

.makeTempFile <- function(name) {
    srcDir <- path.package("AnnotatedMatrix")
    ## When developing with the uncompressed R package, extdata is
    ## still inside inst/ - detect this scenario and accomodate it:
    if (file.exists(file.path(srcDir, "inst")))
        srcDir <- file.path(srcDir, "inst")
    srcDir <- file.path(srcDir, "extdata")
    tmpDir <- tempdir()
    src    <- file.path(srcDir, name)
    if (!file.exists(src)) stop("Failed to identify source file: ", src)
    trg <- NA
    if (dir.exists(src)) {
        ## Source is a directory, set target to the tmp parent
        trg <- tmpDir
        file.copy(src, tmpDir, recursive=TRUE)
    } else {
        ## Basic file
        trg <- file.path(tmpDir, name) 
        file.copy(src, trg)
    }
    sfx    <- gsub('.+\\.', '', name)
    for (sc in sidecarFiles( src )) {
        ## Copy sidecars, too
        file.copy(sc, file.path(tmpDir, basename(sc) ))
    }
    trg
}

#' Extend Data Table
#'
#' Append both new rows and new columns from one data.table to another
#'
#' @details
#'
#' Designed to grow the metadata data.table as both new rows and
#' columns are added, while keeping both rows and columns
#' distinct. There may be a more elegant way to do this with native
#' data.table methods, but I haven't found it, and not for want of
#' trying.
#'
#' Bad things can happen with data.tables if the key column is an
#' integer - when subsetting, these will generally be interpreted as
#' row numbers, which will oftenot correlate to the rows holding the
#' anticipated value, and are also fluid depending on the key(s) being
#' set on the DT (and perhaps on prior operations?) For this reason,
#' the key should be a column of type character.
#'
#' @param x Required, the 'target' data.table
#' @param y Required, the 'source' data.table, which will be added to
#'     or over-write the corresponding rows/columns in x
#' @param key Default 'id', the column to merge on
#'
#' @return A new data.table with new rows and columns merged. Will not
#'     alter input values.
#'
#' @examples
#'
#' library('data.table')
#' ## Again, the key column should be character mode
#' dt1 <- data.table(a = c("1","2","3"), b = c("apple","banana", "cherry"),
#'                   key='a')
#' dt2 <- data.table(a = c("2","4","5"), b = c("BLUEBERRY","DATE","EGGPLANT"),
#'                   c = c("beep","ding","eek"), key='a')
#' extendDataTable(dt1, dt2, key='a')
#'
#' @importFrom data.table rbindlist setkeyv copy
#' @export

extendDataTable <- function (x, y, key='id') {
    ## Find IDs that are already represented
    oldIds <- base::intersect( y[[key]], x[[key]] )
    if (length(oldIds) != 0) {
        ## Add in or update columns for rows that already exist

        ## The update of the rows using rbindlist will generate a
        ## derived object, but altering columns in the manner below
        ## will alter the reference pointed to by x. So we will make a
        ## copy of x to be consistent in behavior (leaving input
        ## untouched):
        x <- data.table::copy(x)
        for (col in colnames(y)) {
            if (col == key) next
            newVals <- y[ oldIds, col, with=FALSE ]
            ## DO WE WANT TO HANDLE NAs SPECIAL?
            ## That is, do we want to NOT overwrite existing values in x
            ## when the y value is NA?
            ## I don't think so ... so fully slotting y into x for these rows
            x[ oldIds, col ] <- newVals
        }
    }
    ## Find any new IDs:
    newIds <- base::setdiff( y[[key]], x[[key]] )
    if (length(newIds) != 0) {
        ## Append the new rows in bulk
        x <- data.table::rbindlist(list(x, y[newIds, , on=key]),  
                                   fill = TRUE, use.names = TRUE)
        data.table::setkeyv(x, key)
    }

    ## The update of the columns (oldIds) should be "in place", but
    ## the rows (newIds) will create a new data.table. Return the result
    x
}
