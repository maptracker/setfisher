
#' Is Defined
#'
#' Broad test returning false for a variety of "empty" values
#'
#' @param x The object to be tested
#' 
#' @aliases is.def
#'
#' @return FALSE for NULL, objects with only NAs, 0x0 matrices and
#'     fields that are classed as 'uninitializedField'. All other
#'     values will be TRUE.
#' 

## Return TRUE for defined objects:
is.def <- function(x) {
    if (is.null(x) || is.empty.field(x)) {
        FALSE
    } else if (is.object(x)) {
        TRUE
    } else {
        ## We will consider matrices "defined" as long as they have
        ## one row OR column. This is to allow certain cached results
        ## to pass after analysis, even if a crossproduct results in
        ## no content. RefClass fields defined as "matrix"
        ## auto-instantiate with 0x0 matrices.
        (is.matrix(x) && (nrow(x) != 0 || ncol(x) != 0)) || any(!is.na(x))
    }
}

## c(0) -> FALSE, c(0,0) -> TRUE
is.something <- function(x) is.def(x) &&
    ( length(x) > 1 || (is.character(x) && x != "") | (is.numeric(x) && x != 0))

## Test for RefClass fields to see if they have been set or not. Unset
## fields are Slot objects of class "uninitializedField"
is.empty.field <- function(x) inherits(x, "uninitializedField")


.parenRE <- function(RegExp, text) {
    ## RegExp in R is not fully fleshed out. This implements a hack
    ## suggested in the internal documentation to allow recovery of
    ## text from multiple parenthetical captures
    m <- lapply(regmatches(text, gregexpr(RegExp, text, ignore.case = TRUE,
                                          perl = TRUE)),
                function(e) {
                    regmatches(e, regexec(RegExp, e, ignore.case = TRUE))
                })
    hitList <- m[[1]]
    ## Return NA if the match failed
    if (length(hitList) == 0) return(NA)
    ## Otherwise, extract out the character vector of the matches, and
    ## exclude the first entry
    hitList[[1]][-1]
    ## Initially written for my first pass at a file registry:
    ## https://github.com/maptracker/GettingAndCleaningData/blob/master/registryManager.R#L153-L166
}

.forwardDeclaration <- function (className = NULL) {
    ## Suggested by Mark Russo 
    if (is.null(className)) return()
    chkParentCls     <- NULL
    try(chkParentCls <- className(className), silent = TRUE)
    if (is.null(chkParentCls)) tmpCls <- setRefClass(className)
}

## I sometimes find it easier to encode lists as string literals
## delimited by newlines
textBlockToVector <- function (x, split="[\n\r]+", trim.white=TRUE) {
    ## Default settings will strip out multiple newlines in a row, but
    ## if you want an empty string you can include a line of just
    ## whitespace (which will default be trimmed away).
    rv <- unlist(strsplit(x, split))
    ## Remove leading and trailing whitespace if requested:
    if (trim.white) rv <- gsub('(^\\s+|\\s+$)', '', rv)
    rv
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
            messsage("File ", n, " requested multiple times")
            next
        }
        cn <- c(cn, n)
        l <- read.table(file, stringsAsFactors=F)
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

