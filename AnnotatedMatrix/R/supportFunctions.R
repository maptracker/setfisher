
## Package-scoped environment to hold items shared across objects
.myPkgEnv <- new.env(parent=emptyenv())
assign(".matrixCache", list(), envir=.myPkgEnv)
## Persistant options internal to annotated matrix
assign(".amParams", list(), envir=.myPkgEnv)


#' Annotated Matrix Parameter
#'
#' Get/Set a parameter associated with your AnnotatedMatrix session
#'
#' @details
#'
#' Manages values that you will likely want to be constant during an
#' analysis, but you may wish to customize between analyses, or
#' between users or organizations.
#'
#' If a requested key is NULL, the function will also check if an R
#' option is set for that key. For example, if you request 'dir', the
#' option 'annotatedmatrixdir' will be checked if the parameter is
#' NULL.
#'
#' @section Recognized Keys:
#'
#' \itemize{
#'
#'   \item \code{dir} - A default directory where matrix files are stored
#'   \item \code{ensemblhost} - Default Ensmebl host domain to use
#' }
#'
#' @param key Required, a character vector of keys to recover. Keys
#'     are case-insensitive, though the names on the returned values
#'     will honor the capitalization of the passed `key` parameter.
#' @param value Optional new value to be assigned to the keys
#' @param default Default \code{NULL}. If the parameter is NULL, and
#'     the corresponding R option is also NULL, this value will be
#'     substitued instead
#' @param as.list Default \code{NULL}. If \code{TRUE}, the return
#'     value will be a named list. \code{FALSE} will return the
#'     unlist()ed values (which may be odd if the parameters have
#'     mixed data types). The default of NULL will return a list if
#'     two or more keys were requested, otherwise the single value
#'
#' @return If as.list=TRUE, a named list. If as.list=FALSE, a single
#'     value if a single key was requested, otherwise an unlist()ed
#'     vector.
#'
#' @importFrom crayon bgYellow
#'
#' @export

annotatedMatrixParameter <- function (key, value=NULL, default=NULL,
                                      as.list=NULL) {
    if (!is.character(key)) stop(crayon::bgYellow(
        "[!] key should be a character, not ", storage.mode(key)))

    lcKey  <- tolower(key)
    params <- get(".amParams", envir=.myPkgEnv)
    if (!is.null(value)) {
        ## User is setting a new value for this key
        for (k in lcKey) {
            params[[ k ]] <- value
        }
        assign(".amParams", params, envir=.myPkgEnv)
    }
    rv <- list()
    for (k in key) {
        v <- params[[ tolower(k) ]]
        if (is.null(v)) {
            ## If we can't find a value set for this key, fallback to
            ## check options() to see if it has been set there, with
            ## an "annotatedmatrix" prefix:
            v <- getOption(paste0("annotatedmatrix",tolower(k) ))
            ## Final fallback to the default value
            if (is.null(v)) v <- default
        }
        rv[[ k ]] <- v
    }
    nkeys <- length(key)
    ## When as.list is not set, it will be TRUE if more than one key
    ## was requested:
    if (is.null(as.list)) as.list <- nkeys > 1
    ## Just return list structure if requested:
    if (as.list) return(rv)
    ## Otherwise, decide if we should use unlist or pluck out a single value:
    if (nkeys == 1) {
        k  <- key[1]
        rv <- rv[[ k ]]
        if (is.vector(rv)) {
            ## Return named
            stats::setNames(rv, k)
        } else {
            ## Naming might be weird, return it 'raw'
            rv
        }
    } else {
        unlist(rv)
    }
}

#' Find Matrices
#'
#' Searches a directory for matrices matching your criteria
#' 
#' @param ns1 Default \code{NULL}. Optional, the "first" namespace
#'     (but will be matched to either the first or second namespace
#'     field)
#' @param ns2 Default \code{NULL}. Options, the "second" namespace,
#'     will be matched to the other field relative to \code{ns1}
#' @param mod Default \code{NULL}. A modifier assigned to some
#'     matrices, often a species name like "Human" or "Mouse".
#' @param type Default \code{NULL}. The type of the matrix, generally
#'     "Map" or "Ontology"
#' @param auth Default \code{NULL}. The authority that provided the
#'     data used to make the matrix, for example "Entrez",
#'     "GeneOntology", "PubMed" etc
#' @param vers Default \code{NULL}. The version of the matrix, eg
#'     "2017-11-05", "GRCh37", etc
#' @param dir Default
#'     \code{getOption("annotatedmatrixdir")}. Required, the folder to
#'     search in on your file system.
#' @param ignore.case Default \code{TRUE}. Should matches ignore case?
#' @param most.recent Default \code{NULL}. If NULL, keep all
#'     results. Otherwise, cluster the matrices by all fields except
#'     version, and keep only the "most recent" version of each
#'     cluster. If \code{most.recent} contains the substring "vers",
#'     this will be done by sorting the Version field of the
#'     data.frame. Otherwise, the Modified field (file system
#'     modification time) will be used.
#' @param recursive Default \code{TRUE}, should subdirectories be
#'     searched as well.
#' @param regexp Default \code{FALSE}. Should matches be performed by
#'     regular expression?
#' @param clear.cache Default \code{FALSE}. File lists are cached for
#'     each \code{dir} requested (as long as \code{recursive} is
#'     TRUE). This allows faster recovery of subsequent searches, but
#'     will prevent discovery of newly-created matrices. Pass a value
#'     of TRUE to assure an explicit search of the directory.
#'
#' @details
#'
#' The expected file format will be one of:
#'
#' \code{<TYPE>@<MOD>-<NS>_to_<NS>@<AUTH>@<VERS>.mtx}
#' \code{<TYPE>@<NS>_to_<NS>@<AUTH>@<VERS>.mtx}
#'
#' ... keeping in mind that the relative order of NS1 and NS2 is
#' irrelevant, and noting that some matricies will not have a modifier
#' (MOD). It will be expected that the suffix be lower case.
#'
#' @return A data.frame with the following fields:
#'
#' \itemize{
#'   \item SubDir - subdirectory relative to \code{dir}
#'   \item Type - The type of conversion matrix
#'   \item Modifier - A modifier, generally a species
#'   \item NS1 - The row namespace / database
#'   \item NS2 - The column namespace / database
#'   \item Authority - The authority for the underlying data
#'   \item Version - The version of the matrix
#'   \item Path - The path to the matrix file, relative to \code{dir}
#'   \item Modified - The file modification time
#' }
#'
#' @seealso \link{annotatedMatrixParameter}
#' 
#' @importFrom CatMisc is.something
#' @export

findMatrices <- function(ns1=NULL, ns2=NULL, mod=NULL, type=NULL, auth=NULL,
                         vers=NULL, dir=annotatedMatrixParameter("dir"),
                         recursive=TRUE, ignore.case=TRUE, most.recent=TRUE,
                         regexp=FALSE, clear.cache=FALSE) {
    if (!CatMisc::is.something(dir)) {
        message("findMatrices(): You must provide the location of your matrix files with the dir parameter")
        return(NA)
    }

    ## I considered building a single, monolithic pattern for
    ## list.files(), but the filters are complex enough that I think
    ## it's best to recover all files, break out the fields, then
    ## filter as a data.frame in code.

    cache    <- get(".matrixCache", envir=.myPkgEnv)
    viaCache <- NULL
    if (recursive && !is.null(cache[[ dir ]]) && !clear.cache) {
        ## Use a previously-cached file list for this directory
        files    <- cache[[ dir ]]
        viaCache <- attr(files, "CacheTime")
    } else {
        ## Find the .mtx files in this folder:
        files <- list.files(path=dir, pattern='\\.mtx$', recursive=recursive)
        if (recursive) {
            ## Cache the result as long as a default recursive search
            ## was run
            attr(files, "CacheTime") <- Sys.time()
            cache[[ dir ]] <- files
            assign(".matrixCache", cache, envir=.myPkgEnv)
        }
    }
    ## Parse out the fields:
    bits  <- CatMisc::parenRegExp("^(.+/)?([^@]+)@(?:(.+?)-)?(.+?)_to_(.+?)@([^@]+)@([^@]+).mtx$", files, unlist=TRUE)
    bits[ bits == "" ] <- NA # Empty string to NA
    ## Structure as matrix, then as data frame
    mat   <- matrix(bits, ncol=7, byrow=TRUE)
    colnames(mat) <- c("SubDir", "Type", "Modifier", "NS1", "NS2",
                       "Authority", "Version")
    df <- as.data.frame(mat, stringsAsFactors=FALSE)
    df$Path <- files

    ## TO DO - What if supplied with multiple values for any of the filters?

    ## Start filtering the df
    usedFilt <- list()
    
    if (!is.null(ns1) || !is.null(ns2)) {
        ## Filter for one or two namespaces
        if (!is.null(ns1) && !is.null(ns2)) {
            ## Testing for two namespaces. ORDER NOT IMPORTANT
            ns1 <- .filePattern(ns1, regexp)
            ns2 <- .filePattern(ns2, regexp)
            df <- df[ (grepl(ns1, df$NS1, ignore.case=ignore.case) &
                       grepl(ns2, df$NS2, ignore.case=ignore.case)) |
                      (grepl(ns1, df$NS2, ignore.case=ignore.case) &
                       grepl(ns2, df$NS1, ignore.case=ignore.case)), ]

            usedFilt$Namespace <- c(ns1, ns2)
        } else {
            ## Testing for just one namespace. Can match to either
            ## rows or columns.
            ns <- ifelse(is.null(ns1), ns2, ns1)
            ns <- .filePattern(ns, regexp)
            df <- df[ grepl(ns, df$NS1, ignore.case=ignore.case) |
                      grepl(ns, df$NS2, ignore.case=ignore.case) , ]

            usedFilt$Namespace <- ns
        }
    }

    if (!is.null(mod)) {
        ## Modifier filtermod
        usedFilt$Modifier <- mod
        mod <- .filePattern(mod, regexp)
        df <- df[ grepl(mod, df$Modifier, ignore.case=ignore.case), ]
    }
    if (!is.null(type)) {
        ## Type filter
        usedFilt$Type <- type
        type <- .filePattern(type, regexp)
        df <- df[ grepl(type, df$Type, ignore.case=ignore.case), ]
    }
    if (!is.null(auth)) {
        ## Authority filter
        usedFilt$Authority <- auth
        auth <- .filePattern(auth, regexp)
        df <- df[ grepl(auth, df$Authority, ignore.case=ignore.case), ]
    }
    if (!is.null(vers)) {
        ## Version filter
        usedFilt$Version <- vers
        vers <- .filePattern(vers, regexp)
        df <- df[ grepl(vers, df$Version, ignore.case=ignore.case), ]        
    }

    ## Add the modified date to the data frame
    df$Modified <- file.info(file.path(dir, df$Path))$mtime
    if (!is.null(most.recent) && !is.na(most.recent[1]) && most.recent[1]) {
        ## Make a synthetic column for aggregation
        df$temp <- paste(df$Type, df$Modifier, df$NS1, df$NS2, df$Authority,
                         sep='@')
        sortCol <- if (grepl("vers", most.recent, ignore.case=TRUE)) {
            message("Most recent version chosen by max() of $Version column.
  If version identifiers do not sort simply, this choice may be incorrect.
")
            "Version"
        } else {
            message("Most recent version chosen by max() of $Modified column.
  If older files were modified on disk recently, you will not recover the most recent data.
")
            "Modified"
        }
        usedFilt$MostRecent <- paste("Via column", sortCol)
        ## Thanks to Scott for reminding me tapply works well here
        
        ## filt will be a named vector; Names will correspond to the
        ## temporary aggregating column, and values will be the "most
        ## recent" value of Version or Modified:
        filt  <- tapply(df[[ sortCol ]], df$temp, max)
        ## Use the named vector as a "lookup hash" to identify rows
        ## that should be kept:
        df    <- df[ df[[ sortCol ]] == filt[ df$temp ], ]

        df$temp <- NULL # Remove the synthetic column
    }
    
    attr(df, "Filter")    <- usedFilt # Note the applied filters
    attr(df, "Directory") <- dir      # Hang the original directory from the df
    attr(df, "CacheTime") <- viaCache # Note if the data were via a cache
    df
}

#' File Pattern
#'
#' Build a query string for findMatrices that includes some wildcards
#'
#' @param qry Required, the query string
#' @param regexp Defaul \code{FALSE}, which will force a match on the
#'     whole string. A value of TRUE will allow substring matches
#' 
#' 
#' @keywords internal

.filePattern <- function (qry, regexp=FALSE) {
    ## Treat underscores and spaces the same:
    qry <- gsub('[_ ]+', '[_ ]+', qry)
    ## Concatenate multiple values into '(foo|bar|bim|...)' format:
    if (length(qry) > 1) qry <- sprintf('(%s)', paste0(qry, collapse='|'))
    ## Match only whole target if requested:
    if (!regexp) qry <- sprintf("^%s$", qry)
    qry
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
#' integer that it can find, and then pick the one that is numerically
#' smallest. The rationale is that "LOC93" was probably annotated
#' earlier (longer ago) than "LOC13200423", and as such is more likely
#' to be the "more common" / "higher expressed" / "better annotated" /
#' "more popular" of the two.
#'
#' @param x A character vector of strings to pick from
#'
#' @return A single string
#'
#' @seealso \link{map}
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
#' @importFrom stats setNames
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
        changes <- stats::setNames(names[diff], goodNames[diff])
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

#' Annotated Matrix Example File
#'
#' Path to a small example file used for demonstrations
#'
#' @details
#'
#' This is just a wrapper around the internal .makeTempFile() function
#' that pulls up a gene symbol -to-> Entrez Gene ID mapping matrix
#' (subset of human genes). The file is factorized and has metadata on
#' both dimensions.
#'
#' @param file Default \code{"Symbol-To-Gene.mtx"}. The file to
#'     recover from inst/
#' @param list Default \code{FALSE}. If TRUE, then list all available
#'     example files
#'
#' @return The path to the file requested, unless list is TRUE, which
#'     will return a listing of all files in \code{extdata/}
#' @examples
#'
#' ef1 <- annotatedMatrixExampleFile()
#' ef2 <- annotatedMatrixExampleFile("Machines.gmt")
#' 
#' @export

annotatedMatrixExampleFile <- function (file="Symbol-To-Gene.mtx", list=FALSE) {
    if (list) {
        srcDir <- file.path(inst.path.package("AnnotatedMatrix"), "extdata")
        rv     <- list.files(srcDir)
        rv[ !grepl('~$', rv) ]
    } else {
        .makeTempFile(file)
    }
}



#' Development Package Path
#'
#' Get the path to package file, including when the package is locally
#' installed from a local installation
#'
#' @details
#'
#' When developing with the uncompressed R package, some files are in
#' different locations. In particular, exdata/ remains inside the
#' inst/ folder. This method detects the presence of "inst" and
#' includes it in the return path
#'
#' @param pkg Default "AnnotatedMatrix", the package name
#'
#' @return The path to the package folder, or the inst/ subfolder if
#'     it exists
#'
#' @seealso \link[base]{path.package}
#'
#' @export

inst.path.package <- function(pkg="AnnotatedMatrix") {
    srcDir <- path.package(pkg)
    ## When developing with the uncompressed R package, extdata is
    ## still inside inst/ - detect this scenario and accomodate it:
    if (file.exists(file.path(srcDir, "inst")))
        srcDir <- file.path(srcDir, "inst")
    srcDir
}

#' Make Temp File
#'
#' Internal method to help manage example files in extdata/
#'
#' @name DOTmakeTempFile
#' @aliases .makeTempFile
#'
#' @details
#'
#' When matrix files are loaded, a .rds version is generated to aid in
#' rapid loading in the future. These files are made alongside the
#' original flat file. This could cause problems if the file is from
#' inst/extdata, so this method takes a requested file and copies it
#' to a temporary location before loading. It will also copy sidecar
#' files used by the 'primary' file.
#'
#' @param name Required, the basename of the file
#' @param pkg Default "AnnotatedMatrix", the package name
#'
#' @seealso \link{sidecarFiles}
#'
#' @examples
#'
#' tmpFile <- AnnotatedMatrix:::.makeTempFile("Symbol-To-Gene.mtx")
#'

.makeTempFile <- function(name, pkg="AnnotatedMatrix") {
    srcDir <- file.path(inst.path.package(pkg), "extdata")
    tmpDir <- tempdir()
    src    <- file.path(srcDir, name)
    if (!file.exists(src)) stop("Failed to identify source file: ", src)
    trg <- NA
    if (dir.exists(src)) {
        ## Source is a directory, set target to the tmp parent
        trg <- file.path(tmpDir, basename(src))
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
            newVals <- y[ oldIds, col, with=FALSE, on=key ]
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

#' Strict Unique
#'
#' Like make.unique, but assure that non-unique names are ALL altered
#'
#' @details
#'
#' code{make.unique} has the property of not altering the first
#' occurance of a non-unique member. That is, the list ('a','b','a')
#' will generate ('a', 'b', 'a.1'). When doing operations that presume
#' or require uniqueness, this is often undesirable. strict.unique
#' will force all non-unique values to be altered.
#'
#' Unlike make.unique, iterative application of the function to the
#' elements vector will NOT create the same output as a single
#' application to the full vector. For its intended use this is not a
#' problem.
#'
#' @param x Required, a character vector
#' @param token Default '#', the string that will go between the
#'     original name and the numeric iterator. So if x had two entries
#'     of "foo", using the default '#' would yields output values of
#'     "foo#1" and "foo#2"
#'
#' @return
#'
#' A character vector where all values are unique, and all non-
#' unique names are altered
#'
#' @seealso \link{make.unique}
#'
#' @examples
#'
#' fowl <-  c("duck", "duck", "goose")
#' 
#' ## The first duck is not modified:
#' make.unique( fowl )
#' 
#' ## All non-unique names will be altered:
#' strict.unique( fowl )
#' 
#' @export

strict.unique <- function (x, token='#') {
    ## Find duplicated values in the input
    dup <- unique( x[ base::duplicated( x ) ] )
    ## Do nothing if there are none
    if (length(dup) == 0) return( x )
    ## Build an iterator for each unque name. It will only be relevant
    ## for the duplicated ones (later in the ifelse() call)
    uNames  <- unique(x)
    counter <- setNames(rep(0L, length(uNames)), uNames)
    ## Make a name vector where every ID is modified with
    ## '#1', '#2' etc
    numName <- vapply( x, function(z)
        sprintf("%s#%d", z, counter[z] <<- counter[z] +1L), "")
    ## Use that modified name only for rows that are not
    ## unique:
    ifelse( is.element( x, dup ), numName, x )
}

#' Printing map() results
#'
#' Print method for mapResult objects generated by AnnotatedMatrix$map()
#'
#' @details
#'
#' $map() normally returns a data.frame with no additional class
#' assigned. If called with format='vector', it will return a named
#' vector with attributes. The attributes end up being printed by
#' print.default, which is generally annoying. This method simply
#' strips attributes before printing.
#' 
#' @param x The object classed as 'mapResult'
#' @param ... Additional parameters, passed on to print.default()
#'
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' locs <- s2e$map(c("MOT8", "GARS"))
#' locs
#'
#' @method print mapResult
#' @export

print.mapResult <- function (x, ...) {
    a <- attributes(x)
    attributes(x) <- list( names = a$names )
    print.default(x, ...)
}

#' Printing Annotated Matrix filter summarys
#'
#' Print method for mapFilter objects generated by
#' AnnotatedMatrix$filterSummary()
#'
#' @details
#'
#' $filterSummary() returns a data frame detailing the number of rows
#' and columns "zeroed-out" by each applied filter. This method
#' formats those results into a more human-friendly format.
#' 
#' @param x The object classed as 'mapFilter'
#' @param color Default NULL, which will read the "useCol" attribute
#'     of x. Logical flagging if output should be colorized with
#'     crayon.
#' @param ... Additional parameters, passed on to message()
#'
#' s2e  <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' s2e$filterSummary( )
#'
#' @importFrom CatMisc is.something
#' @importFrom crayon red yellow magenta cyan
#' @method print mapFilter
#' @export

print.mapFilter <- function (x, color=NULL, ...) {
    nr <- nrow(x)
    uc <- ifelse( is.null(color[1]), any(attr(x, "useCol")), color[1])
    if (nr == 0) {
        ntxt <- "No rows or columns have been filtered out"
        if (uc) ntxt <- crayon::yellow(ntxt)
        message(ntxt)
        return(invisible(0))
    }
    rv   <- 0
    head <- sprintf("%d filter%s removed at least one row or column:",
                    nr,  ifelse(nr == 1, ' has', 's have'))
    if (uc) head <- crayon::red(head)
    txt <- c(head)
    for (i in seq_len(nr)) {
        m   <- x[i, "metric"]
        if (uc) m <- crayon::magenta(m)
        row <- paste(" ", m, ":")
        for (rc in c("Row", "Col")) {
            v <- x[ i, rc ]
            if (!CatMisc::is.something(v)) next
            rv  <- rv + v
            row <- paste(row, v, sprintf("%s%s", rc, ifelse(v == 1, '', 's')))
        }
        txt <- c(txt, row)
        r <- x[ i, "reason" ]
        if (CatMisc::is.something(r)) {
            r <- paste("    #", r)
            if (uc) r <- crayon::yellow(r)
            txt <- c(txt, r)
        }
    }
    sf <- attr(x, "filters")
    if (CatMisc::is.something(sf)) {
        sf <- c("## Applied filters - parsable with $appliedFilters()",
                sprintf("    %s", sf))
        if (uc) sf <- crayon::cyan(sf)
        txt <- c(txt, sf)
    }
    message(paste(txt,  collapse='\n'), ...)
    invisible(rv)
}


#' Small Number Formatter
#'
#' Attempt to make small numbers more human-readable
#'
#' @details
#'
#' Designed for reporting the fraction of non-zero cells in a sparse matrix
#'
#' @param x Required, a number, presumably between 0 and 1
#'
#' @return
#'
#' If x == 0, "none". If x >= 0.001, x as a percentage with two
#' decimal places. Otherwise, x as a fraction with 1 in the numerator,
#' to 4 significant figures.
#'
#' @importFrom CatMisc is.something

smallNumberFormatter <- function(x) {
    x <- x[1]
    ifelse(!CatMisc::is.something(x), "none",
    ifelse(x < 0.001, sprintf("1/%d", signif(1/x, digits=4)),
           sprintf("%.2f%%", 100 * x)))
}

#' Normalize percent
#'
#' Convert text percentages (eg "4%") to numbers based on a denominator
#'
#' @details
#'
#' Used to parse parameters that can be absolutely defined (eg 53) or
#' relatively defined as a percentage (eg "2.5\%"). For the latter will
#' require a denominator to be provided as well
#'
#' @param x Required, a character or numeric vector containing the
#'     values to be evaluated
#' @param denom Default NA. The "denominator" for the calculations,
#'     which will be used when a percentage is provided in x.
#'
#' @importFrom CatMisc parenRegExp
#'
#' @examples
#'
#' ## Numeric strings and percentages can be mixed together; The denom
#' ## parameter only applies to percentages:
#' AnnotatedMatrix:::normalizePercent( c("5", "0.02", "5\%", "0.02%"), 1000)
#' AnnotatedMatrix:::normalizePercent( 1:10 ) # Numbers are 'left alone'
#' 

normalizePercent <- function (x, denom=as.numeric(NA)) {
    isPerc <- CatMisc::parenRegExp("^([0-9\\.]+)%", x)
    ## Can end up with a bunch of coercions to NA inside the ifelse -
    ## this is normal and don't want to bother user with them:
    suppressWarnings( ifelse( is.na(isPerc), as.numeric(x),
                             as.numeric(isPerc) * denom / 100) )
}

#' Current BioMart Version
#'
#' Report the current 'live' version of BioMart
#'
#' @details
#'
#' Gets available Marts and looks for the one named
#' "Ensembl Genes ##", returning the actual value of '##' as the
#' version.
#'
#' @return \code{NA} if there is a problem, otherwise a single integer
#'     value representing the version, with a name representing the
#'     name of the Mart.
#'
#' @param host Will be passed to biomaRt, defaults to the
#'     \link{annotatedMatrixParameter} "ensemblhost", falling back to
#'     www.ensembl.org.
#' @param \dots Passed to \code{listMarts}. In particular, `host` may
#'     be useful to access a preferred Ensembl server.
#'
#' @importFrom biomaRt listMarts
#' @importFrom CatMisc parenRegExp
#' @importFrom stats setNames
#' @importFrom crayon bgYellow
#'
#' @export

biomartVersion <- function (host=annotatedMatrixParameter('ensemblhost', default="www.ensembl.org"), ...) {
    mList     <- biomaRt::listMarts(host=host, ...)
    ## Find the Ensembl Gene version by grepping it out:
    geneVers  <- CatMisc::parenRegExp("Ensembl Genes (\\d+)", mList$version)
    isEnsGene <- !is.na(geneVers)
    rv        <- as.integer(geneVers[isEnsGene])
    ## Should have found only one thing: Throw error if not and return NA
    if (sum(isEnsGene) == 0) {
        message(crayon::bgYellow("Failed to find BioMart Ensembl Gene version"))
        return(NA)
    } else if (sum(isEnsGene) != 1) {
        message(crayon::bgYellow("WEIRD: Multiple BioMart Ensembl Gene versionsions found: ", paste(rv, collapse=', ')))
        return(NA)
    }
    ## Return the version, with mart name as its name
    ##  (should be "ENSEMBL_MART_ENSEMBL")
    stats::setNames(rv, mList$biomart[isEnsGene])
}

#' Current BioMart Genome Versions
#'
#' Return a vector of genome versions in the current version of BioMart
#'
#' @details
#'
#' Recovers the genomic versions used in BioMart, keyed by dataset
#' name associated with each.
#'
#' @param \dots Passed to \link{biomartVersion}
#'
#' @return A vector of genome build tokens (eg "") with dataset names
#'     as vector names
#'
#' @seealso \link{biomartVersion}
#' 
#' @importFrom biomaRt useMart listDatasets
#' @importFrom stats setNames
#' 
#' @export

biomartGenomeVersions <- function (...) {
    vers <- biomartVersion(...)
    if (is.na(vers)) return(NA)
    mart <- biomaRt::useMart( names(vers) )
    x <- biomaRt::listDatasets(mart)
    stats::setNames(x$version, x$dataset)
}

#' BioMart Common Names
#'
#' Dataset-to-Common-Name lookup, providing "human friendly" species names
#'
#' @details
#'
#' Recovers the species common names used in BioMart, keyed by dataset
#' name associated with each.
#'
#' @param \dots Passed to \link{biomartVersion}
#'
#' @return A vector of species common names (eg "Crab-eating macaque")
#'     with dataset names as vector names
#'
#' @importFrom CatMisc parenRegExp
#' @importFrom biomaRt useMart listDatasets
#' @importFrom stats setNames
#' 
#' @export

biomartCommonNames <- function (...) {
    vers <- biomartVersion(...)
    if (is.na(vers)) return(NA)
    mart <- biomaRt::useMart( names(vers) )
    x <- biomaRt::listDatasets(mart)

    ## Parsing eg "Rat genes (Rnor_6.0)"
    ## Defensively anticipating unusual whitespace around what we want
    cNames <- CatMisc::parenRegExp(' *(^.+?) +genes', x$description)

    ## I also want to normalize these names with forms that are used
    ## by NCBI taxonomy common names. Maintining this as a block of
    ## text for easier editting and extension.
    
    ens2ncbiMap <- "
# Ensembl                # Genbank Common name
Rat                      : Norway rat
Fugu                     : Torafugu
Damara mole rat          : Damara mole-rat
Flycatcher               : Collared flycatcher
Kangaroo rat             : Ord's kangaroo rat
Damara mole rat          : Damara mole-rat
Pika                     : American pika
Gorilla                  : Western gorilla
C.savignyi               : Pacific transparent sea squirt
Shrew                    : European shrew
Platyfish                : Southern platyfish
Fruitfly                 : Fruit fly
Lesser hedgehog tenrec   : Small Madagascar hedgehog
Cow                      : Cattle
C.intestinalis           : Vase tunicate
Hyrax                    : Cape rock hyrax
Bushbaby                 : Bush baby
Squirrel                 : Thirteen-lined ground squirrel
Cat                      : Domestic cat
Gibbon                   : Northern white-cheeked gibbon
Naked mole-rat male      : Naked mole-rat
Tarsier                  : Philippine tarsier
Panda                    : Giant panda
Anole lizard             : Green anole
Tree Shrew               : Chinese tree shrew
Chinese softshell turtle : Chinese soft-shelled turtle
Wallaby                  : Tammar wallaby
Guinea Pig               : Domestic guinea pig
Microbat                 : Little brown bat
Coquerels sifaka         : Coquerel's sifaka
"
    for (line in unlist(strsplit(ens2ncbiMap, "\n"))) {
        e2n <- CatMisc::parenRegExp('(.+?)\\s+:\\s+(.+?)\\s*$', line)
        ens <- e2n[1]
        if (is.na(ens) || !is.element(ens, cNames)) next
        cNames[ which(ens == cNames) ] <- e2n[2]
    }
    stats::setNames(cNames, x$dataset)
}
