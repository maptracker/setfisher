#' Parse MatrixMarket File
#'
#' Parses a MatrixMarket file so it can be read by AnnotatedMatrix
#'
#' @details
#'
#' Reads the "normal" matrix data from a MatrixMarket file (the
#' [i,j,x] triples), as well as parses '%' comments to extract
#' metadata usable by AnnotatedMatrix.
#'
#' @param file Required - the path to the .mtx file. Files can be
#'     gzipped, provided they have a .gz suffix.
#'
#' @importFrom CatMisc parenRegExp .flexFilehandle
#' @importFrom ParamSetI selfSplittingString
#' @importFrom Matrix readMM
#' 

## Metadata can be held either in the comments of the MatrixMarket
## file, or as a separate file with a related name. Using a separate
## file is faster (roughly two fold)

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

parse_MatrixMarket_file <- function( file ) {
    rv  <- list()
    ## Begin by simply reading the matrix data
    fhDat <- CatMisc::.flexFilehandle(file)
    if (is.na(fhDat[1])) {
        message("[!]", "MatrixMarket file '", file,"' was not found");
        return(NULL)
    }
    mat <- readMM(fhDat$fh)
    close(fhDat$fh)
    
    ## Now begin parsing comments. Check to see if the file includes
    ## metadata First, find how many rows are comments, then extract
    ## just those rows as a character vector. Doing this in the
    ## (potentially misguided) hope that streaming readLines through
    ## grep and max will minimize RAM overhead for large files
    message("  Reading comment lines...")
    fhDat <- CatMisc::.flexFilehandle(file)
    lastCom <- max(grep("^%", readLines(fhDat$fh)))
    close(fhDat$fh)
    
    fhDat <- CatMisc::.flexFilehandle(file)
    allCom  <- readLines(fhDat$fh, warn = FALSE, n = lastCom)
    close(fhDat$fh)

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
        message("Multiple separators defined in matrix file ",
                file, " : '", paste(sep, collapse="/"), "'")
        return(NULL)
    } else {
        sep <- CatMisc::parenRegExp("Separator\\s+\\'([^\\']+)\\'", allCom[sep])
    }
    fac    <- grep("^LEVELS\\s+\\[.+\\]", allCom)
    if (length(fac) == 0) {
        fac <- NA
    } else if (length(fac) > 1) {
        message("Multiple factor levels defined in matrix file ",
                file, " : '", paste(fac, collapse="/"), "'")
        return(NULL)
    } else {
        facDat    <- CatMisc::parenRegExp("LEVELS\\s+(.+?)$", allCom[fac])
        lvl       <- ParamSetI::selfSplittingString(facDat[1])
        ## Remove leading and trailing whitespace:
        rv$levels <- gsub('^\\s+','', gsub('\\s+$', '', lvl))
    }
    
    ## Get any defaults that have been set
    def <- grep("^DEFAULT\\s+", allCom)
    params <- list()
    for (defl in def) {
        kv <- CatMisc::parenRegExp("DEFAULT\\s+(\\S+)\\s+(\\S.*?)\\s*$",
                                   allCom[defl])
        k <- tolower(kv[1])
        if (is.null( params[[ k ]] )) {
            params[[ k ]] <- kv[2]
        } else {
            message("Ignoring duplicate value for DEFAULT ",k, ' = ', kv[2]);
        }
    }
    rv$params <- params
    
    ## Find the row and column header lines, and optionally Cell metadata
    rcPos <- grep("^(Row|Col|Cell) Name", allCom)
    rcLen <- length(rcPos)
    ## Function to extract each comment row as metadata. strsplit() on
    ## spaces discards trailing spaces, which causes problems for a
    ## space-containing metadata token (eg the default ' :: ') when
    ## the last metadata column is empty. So use CatMisc::parenRegExp
    ## to get the index number instead

    parseRow <- if (is.na(sep)) {
        ## No metadata separator defined
        function (x) CatMisc::parenRegExp("^(\\d+)\\s+(.+)", x)
    } else {
        function (x) {
            indDat <- CatMisc::parenRegExp("^(\\d+)\\s+(.+)", x)
            ## Split the data section on the separator
            c(indDat[1], unlist(base::strsplit(indDat[2], sep)) )
        }
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
        targDat <- CatMisc::parenRegExp("^(\\S+)\\s+(.+?)\\s*$",
                                        allCom[rcPos[i]])
        what    <- targDat[1] # Row / Col
        meta    <- targDat[2]
        if (!is.na(sep))
            meta <- base::strsplit(meta, sep, fixed = TRUE)[[1]]
        meta[1] <- "id" # Generally will be "Name", but normalize to "id"
        mdCol   <- length(meta) - 1;
        hasMD <- if (mdCol > 0) {
            sprintf("and %d metadata column%s", mdCol, ifelse(mdCol==1,"","s"))
        } else {
            "(no metadata)"
        }
        message("  Parsing ", what, " names ", hasMD)
        ## What are the first and last coordinate rows?
        s <- rcPos[i] + 1
        e <- ifelse(i == rcLen, length(allCom), rcPos[i+1] - 1)
        rows <- allCom[s:e]
        ## Remove non-coordinate rows
        isCoord <- grep("^\\d+\\s+\\S", rows)
        ## Make a (possibly) ragged list
        ragged  <- lapply(rows[ isCoord ], parseRow)
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
            ## Reliably building the DT proved non-intuitive. Using
            ## the matrix directly with as.data.table() seems to be
            ## the most reliable mechanism.

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
#'     (rNames) will be taken as names(data). The option to provide
#'     the names independently is included in case some names are not
#'     unique. They will still be made so in the final matrix
#' @param val Default 1. The value to be assigned to non-zero cells in
#'     the matrix.
#' @param file Default NULL. If non-null, will be used to look for
#'     sidecar metadata files.
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

matrixFromLists <- function(data, meta=NULL, listNames=NULL,
                            file=NULL, val=1) {
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
    if (!is.null(meta) && length(meta) > 0) {
        ## Get all the unique IDs for which there is metadata
        ids   <- unique(unlist(lapply(meta, names)))
        ## Make columns
        mList    <- lapply(meta, function(x) x[ids])
        mList$id <- ids
        metadata <- data.table::as.data.table( mList,  key = "id")
        data.table::setkeyv(metadata, "id")
    }
    ## Also load any sidecars:
    if (!is.null(file)) metadata <- parseMetadataSidecar( file, metadata )

    list(matrix = mat, metadata=metadata,
         colChanges = cnDat$changes, rowChanges = cnDat$changes )
}

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

#' Sidecar Files
#'
#' Given a 'main' file, finds any associated / companion files
#'
#' @details
#'
#' Sidecar files are seperate files that carry supporting information
#' for the main data file. In this package they are used to provide
#' additional metadata to file formats that can not store such data
#' internal to the main file
#'
#' @param file Required, the path to the 'main' file
#'
#' @return A character vector of file paths
#'
#' @examples
#' 
#' ## .makeTempFile() is an internal helper function; Don't use normally!
#' lolFile <- AnnotatedMatrix:::.makeTempFile("ListOfLists.inp")
#'
#' sidecarFiles( lolFile )
#'
#' @export

sidecarFiles <- function( file ) {
    fileDir   <- dirname(file)
    ## Remove trailing slashes (eg on directories):
    baseFile  <- gsub("/+$", "", basename(file))
    metafiles <- list.files(fileDir, ignore.case=TRUE,
                            pattern = sprintf("^%s-metadata", baseFile))
    ## Exclude junk (currently text editor stuff)
    metafiles <- metafiles[ ! .isTempFile(metafiles) ]
    ## Make paths absolute to fileDir:
    file.path(fileDir, metafiles)
}

.isTempFile <- function(files) {
    grepl('~$', files) # Emacs temporary files
}


#' Parse Metadata Sidecar
#'
#' Loads metadata from files accompanying the 'main' data
#'
#' @details
#'
#' This package primarily expects matrices to be encoded in
#' MatrixMarket files that include extensive 'in-line' information
#' embedded in comments. In some cases it may be desirable to manage
#' metadata in other files, or you may have file formats other than
#' MTX that can't hold metadata internally.
#'
#' In those cases, one or more 'sidecar' files can be included. This
#' function will expect those files to be named
#' \code{<BASENAME>-metadata<SOMETHING>}, where <BASENAME> is the full
#' path of the 'main' data file, and <SOMETHING> is any other text
#' (including no text).
#'
#' The metadata files should be TSV tabular and include a header
#' row. The function will look for an 'id' column - if not found, it
#' will rename the first column 'id' and use it
#'
#' @param file Required, the path to the primary data file (NOT to the
#'     metadata file)
#' @param metadata Default NULL. Optional data.table that already
#'     contains metadata, presumably from another source. If not
#'     provided, an empty data.table will be made.
#' @param na.strings Default \code{c("NA",'-')}, strings that will be
#'     parsed as NA. A dash is included as it is a common null token
#'     in biological data sets.
#' @param verbose Default TRUE, which will report each found file to
#'     the terminal.
#' 
#' @importFrom crayon white
#' @importFrom data.table data.table fread setkeyv setnames
#'
#' @examples
#'
#' ## .makeTempFile() is an internal helper function; Don't use normally!
#' lolFile <- AnnotatedMatrix:::.makeTempFile("ListOfLists.inp")
#'
#' ## There are two sidecars associated with this file
#' parseMetadataSidecar( lolFile )
#' 
#' 
#' @export

parseMetadataSidecar <- function (file, metadata=NULL,
                                  na.strings=c("NA",'-'), verbose=TRUE) {
    sidecars  <- sidecarFiles( file )
    if (is.null(metadata)) metadata <-
          data.table::data.table( id = character(), key = "id" )
    for (path in sidecars) {
        if (verbose) message("Reading metadata file - ", crayon::white(path))
        scDT <- data.table::fread(path, na.strings=na.strings)
        ## Find the "id" column, or the first one
        scNm    <- names(scDT)
        scIdCol <- which(scNm == "id")
        if (length(scIdCol) == 0) {
            ## No column named 'id'. *Presume* it is the first
            ## one and rename it.
            data.table::setnames(scDT, scNm[1], "id")
        }
        ## Assure that the id column is a character:
        scDT[["id"]] <- as.character(scDT[["id"]])
        data.table::setkeyv(scDT, "id")
        ## Set the key to "id", and then merge into metadata
        metadata <- extendDataTable(metadata, scDT)
        data.table::setkeyv(metadata, "id")
    }
    metadata
}


#' Parse List-of-Lists File
#'
#' Parses a List-of-Lists file so it can be read by AnnotatedMatrix
#'
#' @details
#'
#' This method is somewhat slow; For many files, you may wish to use
#' the Perl script "listOfList_to_MatrixMarket.pl" script in the perl/
#' subdirectory of the AnnotatedMatrix installation.
#'
#' Once the file is parsed, an RDS version will be written "next to"
#' it (same path, plus ".rds"), which will load rapidly, unless a
#' change to the original file triggers an automatic reparsing.
#'
#' @param file Required - the path to the LOL file. Files can be
#'     gzipped, provided they have a .gz suffix.
#'
#' @importFrom CatMisc parenRegExp .flexFilehandle
#' 

parse_ListOfLists_file <- function( file ) {

    data       <- list() # List of character vectors
    nowParsing <- NA     # Name of list currently being parsed
    listLen    <- 0      # Lenght of list currently being parsed
    meta       <- list() # List metadata

    fhDat <- CatMisc::.flexFilehandle(file)
    if (is.na(fhDat[1])) {
        message("[!]", "List-of-list file '", file,"' was not found");
        return(NA)
    }
    while (length(line <- readLines(fhDat$fh, n = 1, warn = FALSE)) > 0) {
        if (line == "") next # Blank line
        if (grepl('^#', line)) {
            ## Comment line
            nameDat <- CatMisc::parenRegExp("^#\\s+LIST\\s+-\\s+(.+?)\\s*$", line)
            nm <- nameDat[1]
            if (!is.na(nm)) {
                ## New list name
                ## Trim up the previous list:
                if (listLen > 0)
                    data[[nowParsing]] <- data[[nowParsing]][ seq_len(listLen) ]
                nowParsing <- nm
                if (is.null( data[[ nowParsing ]] )) {
                    ## First time seen, good
                    data[[ nowParsing ]] <- character(1000)
                    listLen <- 0L
                } else {
                    ## Every LIST entry should be unique
                    stop("List-of-List file: ", file,
                         "\n   has multiple LIST entries with name: ", nm)
                }
            } else {
                ## Check for other metadata, which applies to list only
                metDat <- CatMisc::parenRegExp("^#\\s+([^=]+?)\\s*=\\s*(.+?)\\s*$", line)
                k <- metDat[1]
                if (!is.na(k)) {
                    if (is.null( meta[[ k ]] )) meta[[ k ]] <- character()
                    meta[[ k ]][ nowParsing ] <- metDat[2]
                }
            }
            next # Nothing else to do with comments
        }
        if (grepl('^\\S', line)) {
            ## Presume to be an ID. Extend the character vector if needed:
            if (length( data[[ nowParsing ]] ) == listLen)
                data[[ nowParsing ]] <- c(data[[ nowParsing ]], character(1000))
            listLen <- listLen + 1L
            data[[ nowParsing ]][ listLen ] <- line
        }
    }
    close(fhDat$fh)
    ## Trim up the last list
    if (listLen > 0) data[[nowParsing]] <- data[[nowParsing]][seq_len(listLen)]
    matrixFromLists(data, meta=meta, file=file)
}

#' Parse Text File
#'
#' Parses one or more simple text files into data consumable by AnnotatedMatrix
#'
#' @details
#'
#' Provided one or more files and/or directories, take all files and
#' treat each as a single list
#'
#' Generally you would not call this function directly.
#'
#' @param files Required. One or more file/directory paths (can be a
#'     mixture of directories and files). Files can be gzipped,
#'     provided they have a .gz suffix.
#' @param header Default FALSE. If true, the first row of each file
#'     will be taken as the list name for that list
#' @param rm.suffix Default TRUE, which will remove the file suffix
#'     when the file name is used as a list name
#' @param pattern Default NULL. Optional pattern to use when
#'     extracting files from requested directories
#'
#' @return A list structure usable by AnnotatedMatrix
#'
#' @importFrom CatMisc .flexFilehandle
#' 

parse_Text_file <- function (files, header=FALSE, rm.suffix=TRUE,
                             pattern=NULL) {
    chk   <- files
    files <- character()
    for (file in chk) {
        ## Check each passed argument to see if it is a file or directory
        if (dir.exists(file)) {
            ## The argument is a directory. Grab all files there
            dir   <- file
            dirf  <- list.files(file, ignore.case=TRUE, pattern = NULL)
            ## Exclude junk (currently text editor stuff)
            dirf <- dirf[ ! .isTempFile(dirf) ]
            ## Make paths absolute to the directory:
            files <- c(files, file.path(dir, dirf))
        } else if (file.exists(file)) {
            files <- c(files, file)
        } else {
            message("[?] Requested text file '",file,
                    "' does not appear to exist")
        }
    }
    
    data <- list() # List of character vectors
    for (file in files) {
        fhDat <- CatMisc::.flexFilehandle(file)
        if (is.na(fhDat[1])) {
            message("[!]", "Simple list file '", file,"' was not found");
            next
        }
        
        listname <- if (header) {
            ## Request to take the first line as list name
            readLines(fhDat$fh, n = 1, warn = FALSE)
        } else {
            if (rm.suffix) { fhDat$name } else {basename(file)}
        }
        
        listLen            <- 0L # Lenght of list currently being parsed
        if ( is.null(data[[ listname ]]) ) {
            data[[ listname ]] <- character(1000)
        } else {
            ## We have ended up with the same list name
            ## somehow. Rather than try to figure this out, we will
            ## just extend the existing list
            message("Multiple files have data for '", listname,"'")
            listLen <- length(data[[ listname ]])
        }
        while (length(line <- readLines(fhDat$fh, n = 1, warn=FALSE)) > 0) {
            if (line == "" || grepl("^#", line)) next # Blank line / comment
            if (grepl('\\s', line)) {
                ## Normalize whitespace to single spaces
                line <- gsub('\\s+', ' ', line)
                ## Leading and trailing whitespace removal
                line <- gsub(' $', '', gsub('^ ', '', line))
            }
            if (line != "") {
                if (length( data[[ listname ]] ) == listLen)
                    data[[ listname ]] <- c(data[[ listname ]], character(1000))
                listLen <- listLen + 1L
                data[[ listname ]][ listLen ] <- line
            }
        }
        close(fhDat$fh)
        if (listLen == 0) {
            message("Empty list for '",listname,"' is being ignored")
            data[[listname]] <- NULL
        } else {
            data[[listname]] <- data[[listname]][seq_len(listLen)]
        }
    }
    matrixFromLists(data, file=file)
}

#' Parse GMT File
#'
#' Parses one or more lists from a Gene Matrix Transposed file
#'
#' @details
#'
#' GMT files allow an arbitrary number of lists to be stored, each
#' with a name, description, and one or more members
#'
#' Generally you would not call this function directly.
#'
#' @param file Required, the path to the GMT file
#'
#' @return A list structure usable by AnnotatedMatrix
#'
#' @importFrom CatMisc .flexFilehandle
#' 

parse_GMT_file <- function (file) {
    fhDat <- CatMisc::.flexFilehandle(file)
    names <- character()
    descr <- character()
    data  <- list()
    lnum  <- 0
    while (length(line <- readLines(fhDat$fh, n = 1, warn = FALSE)) > 0) {
        row   <- unlist(strsplit(line, "\t"))
        names <- c(names, row[1])
        descr <- c(descr, row[2])
        lnum  <- lnum + 1
        data[[lnum]] <- row[ -(1:2) ]
    }
    close(fhDat$fh)
    matrixFromLists(data, listNames=names, file=file,
                    meta=list(Description=setNames(descr,names)) )
}

#' RDS is Current
#'
#' Checks if a serialized RDS file is current relative to source files
#'
#' @details
#'
#' For speed, matrix files and their supporting annotation files are
#' serialized to an RDS object. These RDS files will be preferentially
#' used, unless it is observed (by file dates)
#'
#' @param file Required, the path to the 'main' file (not the RDS)
#'
#' @return
#'
#' \itemize{
#' 
#'   \item TRUE - if the RDS object exists and is newer than all
#'         supporting files
#'
#'   \item NA - if the RDS file exists but the original (main) file
#'         does not
#'
#'   \item FALSE - All other situations; No RDS, or RDS is older than
#'         the main file and/or any of the sidecar files
#'
#' }
#' 
#' @export

rdsIsCurrent <- function (file) {
    objFile  <- paste(file,'rds', sep = '.')
    ## If RDS is not present, it's not current either:
    if (!file.exists(objFile)) return( FALSE )
    ## If we don't have the original file, kind of hard to tell:
    if (!file.exists(file)) return( NA )
    ## Original and RDS present. Have any original files changed since
    ## the RDS file was generated? If so, it is no longer current.
    mtime    <- file.mtime(objFile)
    sidecars <- sidecarFiles( file )
    all( mtime >= file.mtime(c(file, sidecars)) )
}
