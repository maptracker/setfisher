## Be quiet when loading package: https://stackoverflow.com/a/8681811

#' SetFisher
#'
#' Enrichment analysis using hypergeometric distribution and adjusting
#' counts for multiple voting
#'
#' @field defQryObj Optional default Query object (SetFisherMatrix)
#' @field defMapObj Optional default Map object (SetFisherMatrix)
#' @field defWorld Optional default World (character vector)
#' @field lastTop Last results from topResults()
#' @field matrixCache Internal storage for SetFisherMatrix objects
#' @field anaCache Internal storage for SetFisherAnalysis objects
#'
#' @importFrom dynamictable dynamictable
#' @importFrom methods setRefClass new
#' @importFrom CatMisc methodHelp
#'
#' @import ParamSetI
#' @import AnnotatedMatrix
#'
#' @export SetFisher
#' @exportClass SetFisher

SetFisher <-
    setRefClass("SetFisher",
                fields = list(
                    defQryObj   = "ANY", ## SetFisherMatrix
                    defMapObj   = "ANY", ## Will cause recursion if set, though
                    defQryWorld = "ANY", ## Default world of IDs
                    matrixCache = "list",
                    lastTop     = "data.frame", # last topResults()
                    anaCache    = "list"
                    ),
                contains = c("ParamSetI")
                )


SetFisher$methods(
    
    ## https://stackoverflow.com/a/13518804
    initialize = function(..., params=NA, matrixCache=list(), help=FALSE ) {
        "Create a new SetFisher object; Invoke with SetFisher(...)"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }
        callSuper(matrixCache=matrixCache, params=params, paramDefinitions="
", ...)
    },

    helpSections = function ( help=FALSE ) {
        "Static list organizing object methods into conceptual sections"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        list(
            "Logging / Messaging" =
                c("message", "dateMessage", "actionMessage", "debugMessage",
                  "err", "verbose", "logText", "showLog"),
            "Utility Methods" = c("tidyTime")
            )
    },

    fieldDescriptions = function(help=FALSE) {
        "A static list of brief descriptions for each field in this object"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        list(
            "log"      = "data.table containing actual logged data",
            "vb"       = "Verbosity flag, set with $verbose()", 
            "EvLogObj" = "Optional 'external' EventLogger object, for log centralization")
    },

    defaultQuery = function (x=NULL, help=FALSE, ... ) {
        "Set/Get an optional default query matrix"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }
        if (!is.null(x)) {
            invisible(defQryObj <<- annotatedMatrix( x, ... ))
        } else if (CatMisc::is.empty.field(defQryObj)) {
            NULL
        } else {
            defQryObj
        }
    },

    defaultMap = function (x=NULL, help=FALSE, ... ) {
        "Set/Get an optional default map matrix"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        if (!is.null(x)) {
            invisible(defMapObj <<- annotatedMatrix( x, ... ))
        } else if (is.empty.field(defMapObj)) {
            NULL
        } else {
            defMapObj
        }
    },

    defaultQueryWorld = function ( x=NULL, help=FALSE ) {
        "Set/Get an optional default query world"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        if (!is.null(x)) {
            invisible(defQryWorld <<- x)
        } else if (is.empty.field(defQryWorld)) {
            NA
        } else {
            defQryWorld
        }
    },

    annotatedMatrix = function ( x=NA, help=FALSE, ... ) {
        "Recover an AnnotatedMatrix, from the internal cache if available"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        if (is.character(x)) {
            ## Presumably a file path
            if (is.null(matrixCache[[ x ]])) {
                ## Not yet in the cache
                matrixCache[[ x ]] <<- AnnotatedMatrix(file=x, ...)
            } else {
                ## Already cached
                matrixCache[[ x ]]
            }
        } else if (inherits(x, "AnnotatedMatrix")) {
            ## Already a matrix object, just return it as-is.
            ## There isn't really a good anchor to
            x
        } else {
            err("Failed to recover matrix - not sure what to do with request of type ", storage.mode(x))
            NULL
        }
    },

    analysis = function (ontology, query=NULL, idmap=NULL,
        name=NULL, queryworld=NULL, help=FALSE, ... ) {
        "Create a SetFisher analysis object or recover one already created"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        numStored <- length(anaCache)
        if (is.numeric(ontology) && ontology > 0 &&
            all.equal(ontology, floor(ontology))) {
            ## If the first argument is an integer number, presume
            ## that an already loaded analysis is being requested by
            ## index.
            ind <- as.integer(ontology)
            if (ind > numStored) {
                err(sprintf("Request for analysis #%d, but only %d available",
                            ind, numStored))
                return(NA)
            }
            return( anaCache[[ ind ]] )
        }
        if (is.character(ontology) && numStored > 0) {
            ## The primary distinguishing aspect of analyses will
            ## generally be their ontologies. See if a passed request
            ## is a file path or name for an already-loaded analysis
            ## ontology.
            pathMatch <- integer()
            nameMatch <- integer()
            for (i in seq_len(length(anaCache))) {
                if (ontology == anaCache[[i]]$ontoObj$file)
                    pathMatch <- c(pathMatch, i)
                if (ontology == anaCache[[i]]$ontoObj$param("name"))
                    nameMatch <-c(nameMatch, i)
            }
            ## If a unique path match was found, return it
            if (length(pathMatch) == 1) return(anaCache[[pathMatch]])
            ## If a unique name match was found, return it
            if (length(nameMatch) == 1) return(anaCache[[nameMatch]])
        }
        qObj <- annotatedMatrix(if (is.null(query)) {
                                    defaultQuery() } else { query })
        if (!CatMisc::is.def(qObj)) err("SetFisherAnalysis objects must define 'query' when created", fatal = TRUE)
        oObj <- annotatedMatrix( ontology )
        if (!CatMisc::is.def(oObj)) err("SetFisherAnalysis objects must define 'ontology' when created", fatal = TRUE)
        mObj <- if (CatMisc::is.def(idmap)) {idmap} else {defaultMap()}
        mObj <- if (CatMisc::is.def(mObj))  {
                    annotatedMatrix(mObj) } else { NA }
        qWld <- if (CatMisc::is.def(queryworld)) {
                    queryworld } else { defaultQueryWorld() }
        
        ## The analysis is the unique combination of the three
        ## matrices. Build a key based on the file paths representing
        ## them, use to assure uniqueness. Note that this means you
        ## can not load multiple analyses from the same files but with
        ## different parameters.
        anaKey <- paste(oObj$file, qObj$file,
                        if (CatMisc::is.def(mObj)) { mObj$file } else {'NA'},
                        collapse = " + ")
        if (is.null(anaCache[[ anaKey ]])) {
            if (is.null(name)) {
                oname <- oObj$param("name")
                if (CatMisc::is.def(oname)) {
                    name <- sprintf("%s Analysis", oname)
                } else {
                    name <- sprintf("Analysis #%d", anaNum)
                }
            }
            ana <- SetFisherAnalysis(
                query = qObj, ontology = oObj, idmap = mObj,
                queryworld=qWld,  setfisher = .self, ...)
            if (!is.null(name)) ana$param("name",name)
            anaCache[[ anaKey ]] <<- ana
        }
        anaCache[[ anaKey ]]
    },

    analysesFromDirectory = function(dir, query=NULL, idmap=NULL,
                                     pattern=NULL, skip=NULL, help=FALSE, ...) {
        "Make a SetFisherAnalysis object for each ontology in a directory"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        files <- list.files(dir, pattern=pattern, full.names=TRUE, ...)
        rv    <- list()
        l     <- 0
        for (f in files) {
            ## If a skip pattern was provided, skip any matching file names
            if (!is.null(skip) && any(grepl(skip, f))) next
            l <- l + 1
            rv[[ l ]] <- NA
            a <- analysis(query=query, ontology=f, idmap=idmap, ...)
            rv[[ l ]] <- a
        }
        rv
    },

    processAll = function ( help=FALSE, ... ) {
        "Filter, process, adjust and report all analyses held by SetFisher"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        for (ana in anaCache) {
            try({
                ana$filter( )
                ana$processListOfLists( )
                ana$adjustResults( )
            })
        }
        topResults( ... )
    },

    analysisNames = function (trim.names=TRUE, help=FALSE, ...) {
        "Takes a set of analysis and shortens their names by removing 'common bits'"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        oNames <- character()
        for (ana in anaCache) {
            aname   <- ana$param("name")
            oNames[ aname ] <- aname
        }
        if (length(oNames) > 1 && trim.names) {
            ## Remove common text from the left and right of each
            ## ontology name. Make certain we don't annhilate a name
            ## doing so.
            prfx <- .sharedSubstring(oNames)
            if (prfx != "") {
                shorter <- gsub(paste('^', prfx, sep=''), '', oNames)
                if (!any(nchar(shorter) == 0)) oNames <- shorter
            }
            sfx <- .sharedSubstring(oNames, left.side = FALSE)
            if (sfx != "") {
                shorter <- gsub(paste(sfx,'$', sep=''), '', oNames)
                if (!any(nchar(shorter) == 0)) oNames <- shorter
            }
        }
        oNames
    },

    topResults = function ( help=FALSE, ...) {
        "Generate a data frame summarizing the top results from all analyses"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

       merged  <- data.frame(term = character(), listName = character(),
                              logPV = numeric())
        
        numOnto <- length(anaCache)
        ## Make note of names of ontologies, to simplify if requested
        oNames <- analysisNames(...)
        for (ana in anaCache) {
            aname   <- ana$param("name")
            capName <- ifelse(numOnto == 1, aname,
                              sprintf("%d ontologies", numOnto))
            tr    <- NULL
            tryCatch({
                tr   <- ana$topResults(add.urls=TRUE, caption.name=capName,... )
                name <- oNames[ ana$param("name") ]
                noHt <- attr(merged, "NoHits")
                if (nrow(tr) == 0) {
                    ## No results
                    ## Do we want to put a placeholder row in the output?
                    attr(merged, "NoHits") <- c(noHt, name)
                } else {
                    ## Set a column with the ontology source
                    tr$Ontology <- name
                    merged <- base::merge(
                        merged, tr, all.y = TRUE, all.x = TRUE)
                    attr(merged, "NoHits") <- noHt
                    attr(merged, "caption") <- attr(tr, "caption")
                }
            }, error = function(e) {
                err(c(paste("Failed to process analysis", aname), e))
                traceback()
            })
            if (is.null(tr)) {
                ##    err(paste("Failed to process analysis", aname))
                ##    traceback()
            }
        }
        lastTop <<- merged[ order(abs(merged$logPV), decreasing=TRUE), ]
        attr(lastTop, "run") <<- TRUE # Success flag
        lastTop
    },

    tsvtable = function(res=NA, file=NA, force=FALSE, help=FALSE, ...) {
        "Convert a top results data frame to a TSV file"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        if (!CatMisc::is.def(res)) {
            ## No results provided
            if (CatMisc::is.something(attr(lastTop, "run"))) {
                res <- lastTop
            } else {
                if (length(anaCache) == 0) {
                    err("Can not generate tsvfile() without at least one analysis")
                    return(NULL)
                }
                res <- processAll( force=force, ... )
            }
        }
        if (!CatMisc::is.def(file)) file <- tempfile("SetFisher-", fileext=".txt" )
        write.table(res, file, quote=FALSE, sep="\t", na="", row.names=FALSE)
        file
    },

    htmltable = function (res=NA, file=NA, caption=NULL, show.attr=FALSE,
        cols=c("term","term Description", "listName", "logPV", "RawLogPV",
            "i","N","i/N","n","W","n/W","Enrich"), ontology.stats=TRUE,
        footer=character(), style=character(), force=FALSE,
        favicon = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4QIHDyEOlqB2hwAAAAh0RVh0Q29tbWVudAD2zJa/AAAA0ElEQVQ4y6WTvQ3CMBSEv1gIUVCkYg0zQpQJaBiEDbKRM0HECLAGjWNEgZBMKHgRkeXYIE666v3p7uwCYLO9aaABat7ogOZyWp/JoJDhI1AGtR6ockuUXC4jtVJq5BbUiXoN4EF7MB6c0HjQAIvU9uXzjjSGEndA5aFSYlgUD7XqchJHD/pIQx8kE5WoxOUKaIGrsP0mgRAaMIATGj7mDTM002ELDAHtAfYebGTYjkkg14YZmlSMI1xigcs9pL+gJp9nDt23CURNJND6c4y5oRc3enezXu9KOQAAAABJRU5ErkJggg==",
        help=FALSE, ...) {
        "Convert a top results data frame to a dynamic HTML table"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }
       
        if (!CatMisc::is.def(res)) {
            ## No results provided
            if (CatMisc::is.something(attr(lastTop, "run"))) {
                res <- lastTop
            } else {
                if (length(anaCache) == 0) {
                    err("Can not generate basichtml() without at least one analysis")
                    return(NULL)
                }
                res <- processAll( force=force, ... )
            }
        }
        if (!CatMisc::is.def(caption)) {
            caption <- attr(res, "caption")
            if (!CatMisc::is.def(caption)) caption <- "Unified report for SetFisher"
        }

        ## Set some default values
        dtParams <- list()
        ## Default column names
        dtParams[[ 'coltitle' ]] <- list(
            term = "The ontology term being tested for enrichment",
            listName = "The name of your submitted list of IDs",
            logPV = "-log10(adjusted hypergeometric p-value), with negative values indicating UNDER-enrichment",
            RawLogPV = "logPV, but before p-value adjustment. Negative values = UNDER-enrichment",
            i = "The number of IDs in your list with the ontology term",
            N = "The total number of IDs in your list",
            n = "The number of IDs in the 'world' with the term",
            W = "The total number of IDs in the 'world'",
            "i/N" =  "The percentage of IDs in your list assigned to the term",
            "n/W" = "The percentage of IDs in the 'world' assigned to the term",
            "Enrich" = "The relative enrichment of the term in your list compared to the world = (i/N)/(n/W)"
            )
        ## Set default gradient information
        pvOpts  <- list(min = -10, max = 10 ) # Basic log p-value gradient
        ## Percent scale from 0-30%
        perOpts <- list(min = 0, max = 30, colors = c("white", "orange"),
                        highstyle = "background-color: red")
        ## Logarithmic scale from 1/10x - 10x
        foldOpts <- list(binVals   = 10 ^ seq(-1,1, length.out=15),
                         highstyle = "color: yellow; background-color: red",
                         lowstyle  = "color: yellow; background-color: blue")
        dtParams[[ 'gradient' ]] <- list(logPV = pvOpts, RawLogPV = pvOpts,
                                       "i/N" = perOpts, "n/W" = perOpts,
                                       Enrich = foldOpts)
        dtParams[[ 'truncate' ]] <- list('*' = 50 )
        dtParams[[ 'factor' ]]   <- list(Ontology = list(r=255, textFilter=TRUE),
                                       listName = list(g=255, textFilter=TRUE))
        ## Colorize term column by the ontology it comes from:
        dtParams[[ 'byFactor' ]] <- list(term = 'Ontology' )
        dtParams[[ 'hide' ]]     <- list(Ontology = TRUE,
                                       'listName Description' = TRUE,
                                       'term Description' = FALSE)
        dtParams[[ 'title'    ]] <- list(listName = 'listName Description',
                                       term = 'Ontology')
        dtParams[[ 'spacemap' ]] <- list(Ontology = '_', term = '_' )
        dtParams[[ 'percent'  ]] <- list("i/N" = TRUE, "n/W" = TRUE)
        dtParams[[ 'fold'     ]] <- list("Enrich" = TRUE)

        ## Add in basic statstics (marble count) columns:
        cnames <- colnames(res)
        iNnW <- match(c("i","N","n","W"), cnames)
        if (all(!is.na(c(iNnW[1], iNnW[2])))) {
            ## User selection statistics
            res[[ "i/N" ]] <- signif(100 * res$i / res$N, 3)
        }
        if (all(!is.na(c(iNnW[3], iNnW[4])))) {
            ## World statistics
            res[[ "n/W" ]] <- signif(100 * res$n / res$W, 3)
            if (all(!is.na(c(iNnW[1], iNnW[2])))) {
                ## Enrichment fold
                res[[ "Enrich" ]] <- signif(res$i * res$W / (res$N * res$n), 3)
            }
        }
        if (ontology.stats) {
            ## WHAT IN THE WORLD WAS I PLANNING HERE?
            
            ## I think these were to be new columns reflecting somehow
            ## on the ontologies? I need to find my notes...
        }
        dynamictable( res, file=file, header=caption, show.attr=show.attr,
                     cols=cols, footer=footer, style=style, favicon=favicon,
                     params=dtParams, ...)
    },

    
    ##     basichtml(sf$lastTop)
    datatable = function( res=NA, caption=NULL, help=TRUE, ... ) {
        "Convert a top results data frame to a data.table"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        if (length(anaCache) == 0) {
            err("Can not generate datatable() without at least one analysis")
            return(NULL)
        }
        if (!CatMisc::is.def(res)) {
            ## No results provided
            if (CatMisc::is.def(lastTop)) {
                res <- lastTop
            } else {
                res <- topResults( ... )
            }
        }
        if (!CatMisc::is.def(caption)) {
            caption <- attr(res, "caption")
            if (!CatMisc::is.def(caption)) caption <- "Unified report for SetFisher"
        }

        ## Arbitrarily take the first analysis object and use it to
        ## call the SetFisherAnalysis$htmltable method
        anaCache[[1]]$datatable( res = res, caption = caption, ... )
    },

    sdrsGrid = function( file=tempfile("SDRS-", fileext=".tsv" ), ... ) {
        ## A TSV file with lists in rows, ontologies in columns, and
        ## logPV values in cells. Somewhat awkward format supporting
        ## internal legacy pipeline.

        ## The file will have three header rows:
        oNames  <- character() # Ontology name, eg 'GO:0006026'
        oDesc   <- character() # Description, eg 'aminoglycan catabolic process'
        oSource <- character() # Ontology data set, eg 'GeneOntology'

        oNames <- analysisNames(...)
        for (ana in anaCache) {
           aname <- ana$param("name")
           res   <- adjustResults( ... )
           
        }
        ## STILL WORKING ON IT
        "Feature not yet implemented"
    },

    show = function ( ) { cat( setFisherText() ) },

    setFisherText = function ( color=TRUE, help=FALSE ) {
        "Generate a compact text summary of the object, used by show()"
        if (help) {
            print( CatMisc::methodHelp(match.call(), class(.self),
                                       names(.refClassDef@contains)) )
            return(NA)
        }

        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        objName <- .self$.selfVarName("mySetFisher")
        colObj  <- doCol(objName, "white")
        msg <- doCol(sprintf("%s base-level object\n",
                                class(.self)[1]), "blue")
        dQ <- if (CatMisc::is.def(defQryObj)) { defQryObj } else { NA }
        dM <- if (CatMisc::is.def(defMapObj)) { defMapObj } else { NA }
        if (CatMisc::is.def(dQ)) {
            name <- doCol(dQ$param("name", default = dQ$file), "yellow")
            msg <- sprintf("%s  %s ( %s$defQryObj ):\n   %s\n", msg,
                           doCol("Default Query", "blue"), colObj,
                           name)
        }
        if (CatMisc::is.def(dM)) {
            msg <- sprintf("%s  Default Map: %s (%s$defMapObj)\n", msg,
                           doCol(dM$param("name", default = dM$file), "red"), colObj)
        }
        numAna <- length(anaCache)
        if (numAna == 0) {
            msg <- sprintf("%s %s\n  Use %s$%s to add analyses\n", msg,
                           doCol("No analyses added yet", "red"), colObj,
                           doCol("analysis( )", "purple"))
        } else {
            msg <- sprintf("%s  %d analys%s in project:\n", msg, numAna,
                           ifelse(numAna == 1, "is", "es"))
            for (i in seq_len(numAna)) {
                ana <- anaCache[[i]]
                msg <- sprintf("%s  %s$analysis(%s)\n", msg, colObj,
                               doCol(i, "red"))
                qObj <- ana$queryObj
                if (!identical(qObj, dQ)) msg <-
                    sprintf("%s    Query: %s\n", msg, doCol(
                        qObj$param("name", default = "-No name-"), "blue"))
                if (CatMisc::is.def(ana$mapObj) && !identical(ana$mapObj, dM)) msg <-
                    sprintf("%s      Map: %s\n", msg, doCol(
                        ana$mapObj$param("name", default = "-No name-"), "red"))
                msg <- paste(msg, ana$.ontologySummaryText(pad = "  ", color=color))
            }
        }


        msg <- paste(c(msg, doCol("Object methods:\n", "red")))
        msg <- paste(c(msg,sprintf("  %s$%s : %s\n", colObj, doCol("matrix( file )", "purple"), "Load a matrix file")))
        msg <- paste(c(msg,sprintf("  %s$%s : %s\n", colObj, doCol("analysis( ontolog, query, idmap)", "purple"), "Define an analysis")))
        msg <- paste(c(msg,sprintf("  %s$%s : %s\n", colObj, doCol("log", "purple"), "Action log")))

        msg
    })


## This was more elegant than my vapply solution
## From ?strsplit, as pointed out by @josh-obrien :
##   https://stackoverflow.com/a/13613183
.strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev),
                                  paste, collapse="")

#' Shared Substring
#'
#' Find a common substring on the left or right of a set of strings
#' 
#' @name dotSharedSubstring
#'
#' @details
#'
#' Will identify characters common to all provided strings on either
#' the left or right side. Used for trimming away names to produce
#' more compact reports.
#'
#' @param x Required, a character vector of strings
#' @param left.side Default \code{TRUE}, which will process the left
#'     side of the strings. If FALSE it will ... wait for it
#'     ... process the \strong{RIGHT} side of the string.
#'
#' @return A string representing the common part of all strings.
#'
#' @keywords internal

.sharedSubstring <- function (x, left.side=TRUE) {
    ## Lightly adapted from Bioconductor source code for lcPrefix:
    ## https://github.com/Bioconductor-mirror/Biobase/blob/master/R/strings.R
    ##    License: Artistic-2.0
    ## ... as suggested by @rpierce : https://stackoverflow.com/a/33737045
    
    nc <- nchar(x, type="char")
    ## If we are looking at the right side (suffix), reverse all strings:
    if (!left.side) x <- .strReverse(x)
    
    share <- ""
    for (i in 1:min(nc)) {
        ss <- substr(x, 1, i)
        if (any(ss != ss[1])) {
            share <- substr(x[1], 1, i-1)
            break
        }
    }
    if (!left.side) share <- .strReverse(share)
    unname(share)
}
