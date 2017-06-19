#' SetFisher Matrix
#'
#' Annotated sparse matrix for capturing query lists, identifier
#' mappings and ontology lookups
#'
#' @field setfisher Pointer back to the parent SetFisher object
#' @field file Path to file the matrix was loaded from
#' @field fromRDS Logical, true if the loaded file was an RDS object
#' @field log The SetFisherLogger object holding log (activity)
#'     entries
#' @field matrixRaw The Matrix object as loaded from the file
#' @field matrixUse The Matrix object after any filters are applied
#' @field matrixMD data.table holding metadata associated with the
#'     matrix
#' @field levels Character array of level names for factor matrices
#' @field rowChanges Named character vector of any row names that
#'     needed changing. Values are the original name, names are the
#'     names after processing with make.names() (if valid = TRUE) or
#'     make.unique() if (valid = FALSE)
#' @field colChanges As per rowChanges, but for column names
#'
#' @import Matrix
#'
#' @export SetFisherMatrix
#' @exportClass SetFisherMatrix
#' 
#' @include SetFisherLoggerI.R
#' @include SetFisherParamI.R
#' @include SetFisherUtils.R
#'

SetFisherMatrix <-
    setRefClass("SetFisherMatrix",
                fields = list(
                    setfisher  = "SetFisher",
                    file       = "character",
                    fromRDS    = "logical",
                    log        = "SetFisherLogger",
                    matrixRaw  = "dgTMatrix",
                    matrixUse  = "ANY", # dgTMatrix
                    matrixMD   = "data.table",
                    levels     = "character",

                    ## If row or column names need to be remapped:
                    rowChanges = "character",
                    colChanges = "character"
                    ),
                contains = c("SetFisherLoggerI", "SetFisherParamI")
                )


SetFisherMatrix$methods(
    
    initialize = function(... , setfisher = NA, file = NA, param = NA ) {
        if (!is.def(setfisher)) return()
        if (!is.def(setfisher)) stop("SetFisherMatrix entries should be created from a SetFisher object. Provided = ", setfisher)
        log <<- setfisher$log # Set here so messaging works
        if (!is.def(file)) err("SetFisherMatrix must define 'file' when created", fatal = TRUE)
        callSuper(..., setfisher = setfisher, file = file)
        fromRDS   <<- FALSE
        .self$.setParamDefs("
Name        [character] Optional name assigned to the matrix
Description [character] Optional description for the matrix
RowDim      [character] Optional name for the row dimension
ColDim      [character] Optional name for the column dimension
RowUrl [character] Optional base URL for row names (%s placeholder for name)
ColUrl [character] Optional base URL for column names (%s placeholder for name)
")
        .self$.setParamList( param )
        .readMatrix( ... )
        matrixUse <<- NULL
    },
    
    matrix = function( raw = FALSE) {
        if (raw | !is.def(matrixUse)) { matrixRaw } else { matrixUse }
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
        if (!is.def(file)) err("SetFisherMatrix objects must define 'file' when created", fatal = TRUE)
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
            if (!file.exists(file)) err(c("Can not make SetFisherMatrix",
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
            } else {
                err(c("Can not make SetFisherMatrix - unrecognized file type: ",
                      file), fatal = TRUE)
            }
            ## Stub DT if no metdata was defined:
            if (!is.something(rv$metadata))
                rv$metadata <- data.table( id = character(), key = "id" )
            dateMessage(paste("Serializing matrix to file",
                              colorize(objFile,"white")), prefix = "  ")
            saveRDS(rv, objFile)
        }
        matrixRaw  <<- rv$matrix
        matrixMD   <<- rv$metadata
        if (!is.null(rv$rowChanges)) rowChanges <<- rv$rowChanges
        if (!is.null(rv$colChanges)) colChanges <<- rv$colChanges
        if (is.def(rv$levels)) levels <<- rv$levels
        ## Set default parameters, without clobbering any already set
        if (is.def(rv$params)) .setParamList(rv$params, clobber = FALSE)

        ## Numeric conversion to prevent integer overflow on product
        cnum       <- as.numeric(ncol(matrixRaw))
        rnum       <- as.numeric(nrow(matrixRaw))
        nnz        <- nnzero(matrixRaw)
        actionMessage(sprintf("%d x %d matrix, %.2f%% non-zero",
                              rnum, cnum, 100 * nnz / (cnum * rnum)),
                      prefix = "  ")
        rv
    },

    .unique.names = function( names = character(), rpt = NULL,
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
                msg <- colorize(paste(num,rpt,"required alteration"),
                                bgcolor = "cyan")
                if (num <= 20) msg <- c(msg, vapply(seq_len(length(changes)),
                        function (i) {
                            sprintf("  '%s' -> '%s'",changes[i],
                                    names(changes)[i])
                        }, ""))
                message(msg, collapse = "\n")
            }
        }
        list( names = goodNames, changes = changes )
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
                nameDat <- .parenRE("^#\\s+LIST\\s+-\\s+(.+?)\\s*$", line)
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
        rnDat         <- .unique.names(names(inames[1:icnt]),
                                             "Row Names")
        cnDat         <- .unique.names(names(jnames[1:jcnt]),
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
        } else if (!is.something(listname)) {
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
        
        rnDat         <- .unique.names(names, "Row Names")
        cnDat         <- .unique.names(listname, "List Name")
        rownames(mat) <- rnDat$names
        colnames(mat) <- cnDat$names
        list( matrix = mat,
             colChanges = cnDat$changes, rowChanges = cnDat$changes )
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
            sep <- .parenRE("Separator\\s+\\'([^\\']+)\\'", allCom[sep])
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
            facDat <- .parenRE("LEVELS\\s+\\[([^]]+)\\]\\[([^]]+)\\]",
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
            kv <- .parenRE("DEFAULT\\s+(\\S+)\\s+(\\S.*?)\\s*$", allCom[defl])
            params[[ tolower(kv[1]) ]] <- kv[2]
        }
        rv$params <- params
        parseRow <- function (x, sep = NA) {
            ## strsplit() on spaces discards trailing spaces, which
            ## causes problems for a space-containing metadata token
            ## (eg the default ' :: ')when the last metadata column is
            ## empty. So use .parenRE to get the index number instead:

            indDat <- .parenRE("^(\\d+)\\s+(.+)", x)
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
        metadata <- data.table( id = character(), key = "id" )
        dimNames <- list( Row = NULL, Col = NULL )
        dimChngs <- list( Row = NULL, Col = NULL )
        for (i in seq_len(rcLen)) {
            ## Is this specifying Row or Col, and what are the headers?
            targDat <- .parenRE("^(\\S+)\\s+(.+?)\\s*$", allCom[rcPos[i]])
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
            colnames(cmat) <- .unique.names(
                c(tmpIndName, meta), paste("Metadata headers for",what))$names
            rownames(cmat) <- .unique.names(cmat[, "id" ])$names
            ## The first column are indices in the sparse matrix
            inds    <- as.integer(cmat[, 1])
            ## The second column holds the names
            names <- character()
            names[ inds ] <- cmat[, 2]
            mnn <- .unique.names(names,paste(what, "Names"))
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
                tmp <- as.data.table( cmat[, -1, drop = FALSE], key = "id")
                rownames(tmp) <- .unique.names(cmat[, "id" ])$names
                ## I never did find a way to add *ROWS* by reference in DTs...
                metadata <- rbindlist(list(metadata, tmp), fill = TRUE)
                setkey(metadata, "id") ## Merged DTs don't carry over key
            }
        }
        ## Set the dimension names
        rdn <- params$rowdim
        names(dimNames)[1] <- ifelse(is.something(rdn), rdn, "")
        cdn <- params$coldim
        names(dimNames)[2] <- ifelse(is.something(cdn), cdn, "")
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
            scDT <- fread(mf)
            ## Find the "id" column, or the first one
            scNm <- names(scDT)
            scIdCol <- which(scNm == "id")
            if (length(scIdCol) == 0) {
                ## No column named 'id'. *Presume* it is the first
                ## one and rename it.
                setnames(scDT, scNm[1], "id")
            }
            setkey(scDT, "id")
            ## Set the key to "id", and then merge into metadata
            metadata <- rbindlist(list(metadata, scDT), fill = TRUE,
                                  use.names = TRUE)
            setkey(metadata, "id")
        }
        rv$matrix     <- mat
        rv$metadata   <- metadata
        rv$colChanges <- dimChngs$Col
        rv$rowChanges <- dimChngs$Row
        rv
     },
    
    metadata = function ( id = NULL, key = NULL) {
        if (!is.something(key)) {
            ## No column specified
            if (!is.something(id)) return( NULL )
            ## Return a subset of the data table:
            rv <- matrixMD[ id, ]
            setkey(rv, "id")
            return( rv )
        } else if (!is.something(id)) {
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

    show = function (...) { cat( .self$.showText(...) ) },

    .showText = function ( pad = "", useObj = NULL, fallbackVar = NULL,
                          compact = FALSE, color=TRUE ) {
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
        if (is.something(levels) && !compact) {
            ## Report the factor levels
            indent <- "\n      ";
            lvl <- paste(doCol(strwrap(paste(levels, collapse = ', '),
                                          width = 0.7 * getOption("width")),
                                  "purple"), collapse = indent)
            msg <- sprintf("%s    %s%s%s\n", msg, doCol(
                "Values are factors:", "yellow"), indent, lvl)
        }
        mNm <- attr(matrixRaw, "matrixName")
        name <- param("name")
        if (is.something(name))
            msg <- sprintf("%s    Name: \"%s\"\n", msg, doCol(name, "white"))
        desc <- param("description")
        if (is.something(desc) && !compact)
            msg <- sprintf("%s    Desc: %s\n", msg, doCol(desc, "white"))
        
        ## If an external "used" matrix is not provided, use internal field:
        if (!is.def(useObj)) useObj <- matrixUse
        dimNames <- names(dimnames(matrixRaw))
        mat <- if (is.null(useObj)) { matrixRaw } else { useObj }
        nr  <- nrow(mat)
        nc  <- ncol(mat)
        nz  <- nnzero(mat)
        msg <- sprintf("%s  %8d %s", msg, nr, ifelse(is.something(
            dimNames[1]),dimNames[1],"rows"))
        rN  <- rownames(mat)
        if (is.def(rN))
            msg <- sprintf("%s eg: %s", msg, doCol(substr(paste(rN[1:pmin(3,length(rN))], collapse = ', '), 1, 60), "cyan"))
        
        msg <- sprintf("%s\n  %8d %s", msg, nc, ifelse(is.something(
            dimNames[2]),dimNames[2],"cols"))
        cN  <- colnames(mat)
        if (is.def(cN))
        msg <- sprintf("%s eg: %s", msg, doCol(substr(paste(cN[1:pmin(3,length(cN))], collapse = ', '), 1, 60), "cyan"))

        perPop <- sprintf("%.3g%%", 100 * nz / (nr * nc))
        msg <- sprintf("%s\n  %d non-zero cells (%s)", msg, nz,
                       doCol(perPop, "red"))
        msg <- sprintf("%s\n", msg)
        ## Note if we had to remap one or more R/C names
        changed <- character()
        if (is.def(rowChanges)) changed <- c(changed, doCol(
            paste(length(rowChanges), "Rows ($rowChanges)"), "yellow"))
        if (is.def(colChanges)) changed <- c(changed, doCol(
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
            msg   <- paste(c(msg, doCol(sprintf("  %d IDs have metadata assigned in up to %d keys:\n", mdRow, length(mcols)), "blue"), collapse = ""))
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
        if (is.def(pad)) msg <-
            paste(c(sprintf("%s%s", pad, unlist(base::strsplit(msg, "\n"))),
                    "\n"), collapse = "\n")
        ## Ending up with extraneous newlines at end
        sub("\n+$", "\n", msg)
    }
    
)
