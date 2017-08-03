
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
#' @export

annotatedMatrixExampleFile <- function () {
    .makeTempFile("Symbol-To-Gene.mtx")
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
#' @return The path to the package folder, or the inst/ subfolder if
#'     it exists
#'
#' @seealso \link{path.package}
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
#' tmpFile <- AnnotatedMatrix::.makeTempFile("Symbol-To-Gene.mtx")
#'

.makeTempFile <- function(name, pkg="AnnotatedMatrix") {
    srcDir <- file.path(inst.path.package(pkg), "extdata")
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
        sprintf("%s#%d", z, counter[z] <<- counter[z] +1), "")
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
#' s2e$filterSummary(fs)
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
        sf <- c("## Applied filters", sprintf("    %s", sf))
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
