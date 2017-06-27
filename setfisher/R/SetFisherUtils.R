
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

