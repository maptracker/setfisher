
logEtolog10 <- log(10) # Natural to base-10 conversion factor

#' SetFisher Analysis
#'
#' Class object associating the components of an enrichment analysis
#' (query lists, ontology, and optional mapping matrix) with the
#' results
#'
#' @field setfisher Pointer back to the parent SetFisher object
#' @field queryObj AnnotatedMatrix representing the query list(s)
#' @field queryWorld Character verctor describing all objects in the
#'     query namespace. Will default to all unique rows in the mapping
#'     matrix, or all unique rows in the query matrix if no mapping
#'     matrix is present.
#' @field activeQueries Character vector of normalized IDs
#'     representing the world. Will either be taken from
#'     \code{queryWorld}, or dynamically calculated from filtered /
#'     intersected matrices.
#' @field ontoObj AnnotatedMatrix representing the mapping of IDs
#'     (from the mapping matrix if present, otherwise from the query
#'     matrix) to ontology terms
#' @field mapObj Optional AnnotatedMatrix mapping IDs in the Query
#'     matrix to those used by the Ontology matrix
#' @field mapWeights Weight matrix used to accomodate multiple voting
#'     from the Query namespace to the Ontology namespace. Taken from
#'     the filtered (non-raw) Mapping Matrix if defined, or as a
#'     simple identity matrix from the Query Matrix if not.
#' @field ontoUse Boolean matrix indicating allowed associations
#'     between IDs (either input IDs, or mapped IDs if a Mapping
#'     matrix is used) and ontology terms. The matrix is filtered to
#'     align with mapWeights.
#' @field queryOnto A direct Query-to-Ontology matrix that takes into
#'     account the existence or absence of a Mapping Matrix.
#' @field resultRaw 3D array holding raw p-values from phyper()
#'     calculations, as well as i and n values for each calculation.
#' @field resultAdj p.adjust() processed values from resultRaw
#' @field lastTop data.frame of the most recent values from
#'     topResults()
#' @field logThresh log10 value of the user's p-value significance
#'     threshold
#' @field roundUp Integer representing the denominator of fractional
#'     counts that should be rounded up to a full count. That is, a
#'     value of 4 indicates that any fractional count >= 1/4 should
#'     ultimately be considered a full count.
#' @field pseudoRound This is a pseudo-count based off of roundUp. It
#'     is added to values before they are round()ed.
#' @field ontoNames Vector of all valid ontology terms following
#'     filtering
#' @field ontoSize Total number of Ontology terms surviving filtering
#'     = length(ontoNames).
#' @field okInput Standardized Query names, following matrix pruning
#' @field idCount Fractional totals for each ID used in the Ontology
#'     matrix, after pruning
#' @field ontoCount Fractional total ID counts for each Ontology term
#'     following pruning
#' @field worldSize Integer count of total number of IDs in the world,
#'     following pruning
#' @field bogusList Used for debuging - holds the most recent list to
#'     fail completely when run through processList()
#' @field workSpace Internal temporary structure used by topResult()
#'
#' @importFrom dynamictable dynamictable
#' @importFrom methods setRefClass new
#'
#' @import AnnotatedMatrix
#' @import ParamSetI
#' 
#' @export SetFisherAnalysis
#' @exportClass SetFisherAnalysis
#'



### Some fields
## ontoSize  Number of ontology terms that survived trimming
## okInput   Vector of trimmed query terms that can be used for analysis
## idCount   Vector of fractional representations of each (mapped) gene
##    was geneCount
## worldSize Rounded sum of geneCount = m+n (W) the size of the full set
## ontoCount Rounded integer vector of genes assigned to each ontology
## ontoNames Vector of names for all surviving ontology terms


SetFisherAnalysis <-
    setRefClass("SetFisherAnalysis",
                fields = list(
                    setfisher     = "SetFisher",
                    queryObj      = "ANY", # "AnnotatedMatrix",
                    queryWorld    = "character",
                    activeQueries = "character",

                    discardedIDs  = "list",

                    ontoObj       = "ANY", # "AnnotatedMatrix",
                    mapObj        = "ANY", # AnnotatedMatrix
                    mapWeights    = "dgTMatrix",
                    modState      = "numeric",

                    ontoUse       = "ANY", # dgTMatrix
                    queryOnto     = "dgTMatrix", # Direct Qry -> Onto
                    
                    resultRaw     = "array", # [list x ontology x HGD]
                    resultAdj     = "array", # adjusted results
                    lastTop       = "data.frame", # last topResults()
                    
                    logThresh     = "numeric", # -log10 value of threshold
                    
                    ## Fractions >= 1/roundUp will be rounded up. The
                    ## default of 2 means that an input query (eg an
                    ## Affy probe set) that maps to two target IDs (eg
                    ## an Entrez gene) can "count" if only 1 out of 2
                    ## are positive.
                    roundUp        = "integer",
                    ## pseudoRound is just (roundUp - 1) / roundUp
                    pseudoRound    = "numeric",
                    ## qualities of the analysis or world as a whole:
                    ontoSize       = "integer", # Num of ontology terms
                    okInput        = "character", # 'legal' input IDs
                    ontoNames      = "character", # Names of all filtered ontos
                    idCount        = "numeric", # Fractional count for targ IDs
                    ontoCount      = "integer", # Total IDs for each onto
                    worldSize      = "integer", # Total target IDs in world
                    outputDim      = "integer", # Used to 'align' matrices
                    inOutFunc      = "list",    # Simple row/col -> in/out funcs

                    ## Failed list - used for debugging
                    bogusList      = "character",
                    ## scratch variable, to suppress complaints about
                    ## '<<-' needed in an apply() loop
                    workSpace      = "integer"
                    
                    ),
                contains = c("ParamSetI")
                )

SetFisherAnalysis$methods(
    
    initialize = function(setfisher = NULL,
        ontology = NULL, query = NULL, idmap = NULL, queryworld = NA,
        round = 2L, params = NA, ... ) {
        "RefClass initialization method to make new SetFisherAnalysis Object"
        
        usingMethods("param")
        if (!CatMisc::is.def(setfisher)) stop("SetFisherAnalysis entries should be created from a SetFisher object")
        setfisher <<- setfisher
        modState  <<- c(-1, -1, -1) # Modified State: Qry, Map, Onto
        log <<- setfisher$log # Set here so messaging works. MAY NOT BE NEEDED
        
        if (!CatMisc::is.def(query)) log$err(
            "SetFisherAnalysis must define 'query' when created", fatal=TRUE)
        queryObj <<- query
        
        if (!CatMisc::is.def(ontology)) log$err(
            "SetFisherAnalysis must define 'ontology' when created", fatal=TRUE)
        ontoObj <<- ontology
        
        queryWorld <<- if (CatMisc::is.def(queryworld)) {
            if (length(queryworld) == 0) {
                err("Query world was set to an empty vector - ignoring")
                as.character(NA)
            } else {
                .standardizeId(queryworld)
            }
        } else {
            as.character(NA)
        }
        
        roundLevel( round )
        ## Set some parameters
        mapObj <<- if (CatMisc::is.def(idmap)) {
            ## An optional mapping matrix has been provided
            idmap
        } else {
            NA
        }
        ## Initialize structure to hold IDs that fall out due to
        ## failure to map, restricted worlds, or lack of information:
        nec <- stats::setNames(character(), character()) # Empty named char vec
        discardedIDs <<- list(queryID=nec,
                              ontologyID=nec,
                              ontologyTerm=nec )
        threshold(0.05)
        callSuper(..., params=params, paramDefinitions="
Name         [character] Optional name for the analysis
" )
        ##print(str(query))
        ##print(str(queryObj))
        .alignMatrices()
    },
    
    .alignMatrices = function (raw=FALSE) {
        "Internal method that identifies dimension overlap between the matrices"
        ## We're going to just do this at initialization. This is a
        ## fairly fundamental relationship between the matrices, and
        ## if it is sensitive to later operations there's something
        ## fundamentally wrong/messy with the matrix IDs
        
        outputDim  <<- as.integer(c(NA, NA, NA))

        ## A variety of ways this alignment can fail. Try to collect
        ## informative feedback for user if we can't line everything
        ## up.
        issues <- character()
        if (queryObj$nnZero() == 0) issues <- c(issues, "Empty Query matrix.")
        if (ontoObj$nnZero() == 0) issues <- c(issues, "Empty Ontology matrix.")
        if (CatMisc::is.def(mapObj)) {
            ## Three matrices to align
            if (mapObj$nnZero() == 0) issues <- c(issues, "Empty Map matrix.")
            sdQM   <- queryObj$sharedDimensions( mapObj, raw=raw )
            sdMO   <- mapObj$sharedDimensions( ontoObj, raw=raw )
            if (is.na(sdQM[1])) issues <- c(issues, "No overlap found between Query and Mapping.")
            if (is.na(sdMO[1])) issues <- c(issues, "No overlap found between Mapping and Ontology.")
            if (length(issues) == 0 && sdQM[2] == sdMO[1]) issues <- c(issues, "Mapping matrix is not bridging Query and Ontology - it looks like it might not be needed?")
            if (length(issues) > 0) {
                message(c("Failed to align Query, Map and Ontology:", issues),
                        prefix = "[WARN]", bgcolor = "yellow", color = "blue")
            } else {
                outputDim <<- c(ifelse(sdQM[1] == 1L, 1L, 2L),
                                sdMO[1],
                                ifelse(sdMO[2] == 1L, 2L, 1L) )
            }
        } else {
            ## Just a set of queries and an ontology
            sdQO <- queryObj$sharedDimensions( ontoObj, raw=raw )
            if (is.na(sdQO[1])) issues <- c(issues, "Failed to find common ids between Query and Ontology")
            if (length(issues) > 0) {
                message(c("Failed to align Query with Ontology:", issues),
                        prefix = "[WARN]", bgcolor = "yellow", color = "blue")
            } else {
                outputDim <<- c(ifelse(sdQO[1] == 1L, 1L, 2L),
                                as.integer(NA),
                                ifelse(sdQO[2] == 1L, 2L, 1L) )
            }
        }

### Set up some callbacks to streamline input/output requests. The
### matrices will have an inherent "input to output", "left to right"
### polarity that will be defined by the "namespace" overlap
### determined above. The package will be operating presuming queries
### "entering" ("input") from the Query matrix, and enrichment results
### "exiting" ("output") from the Ontology matrix, possibly transiting
### through a Mapping matrix. We are not enforcing that input needs to
### be rows (or columns) so the little functions below transparently
### wrap the choice or row or column after the above alignment has
### occurred.

        ## "Null" functions that just return a single typed-NA when
        ## the dimension could not be determined:
        nullChrFunc <- function(...) as.character(NA)
        nullIntFunc <- function(...) as.integer(NA)
        nullLogFunc <- function(...) as.logical(NA)
        
        ## q = query, m = map, o = ontology
        ## i = input, o = output
        ## n = names, c = counts, p = populated
        inOutFunc  <<- list(qin=nullChrFunc,
                            qon=nullChrFunc,
                            qic=nullIntFunc,
                            qoc=nullIntFunc,
                            qip=nullLogFunc,
                            qop=nullLogFunc,
                            
                            min=nullChrFunc,
                            mon=nullChrFunc,
                            mic=nullIntFunc,
                            moc=nullIntFunc,
                            mip=nullLogFunc,
                            mop=nullLogFunc,
                            
                            oin=nullChrFunc,
                            oon=nullChrFunc,
                            oic=nullIntFunc,
                            ooc=nullIntFunc,
                            oip=nullLogFunc,
                            oop=nullLogFunc)
        
        if (!is.na(outputDim[1])) {
            ## We can make some callbacks for the query matrix
            if (outputDim[1] == 1L) {
                ## Output is represented by rows
                inOutFunc$qin <<- function(...) queryObj$cNames(...)
                inOutFunc$qon <<- function(...) queryObj$rNames(...)
                inOutFunc$qic <<- function(...) queryObj$cCounts(...)
                inOutFunc$qoc <<- function(...) queryObj$rCounts(...)
                inOutFunc$qip <<- function(...) queryObj$populatedColumns(...)
                inOutFunc$qop <<- function(...) queryObj$populatedRows(...)
            } else {
                ## Output is represented by columns
                inOutFunc$qin <<- function(...) queryObj$rNames(...)
                inOutFunc$qon <<- function(...) queryObj$cNames(...)
                inOutFunc$qic <<- function(...) queryObj$rCounts(...)
                inOutFunc$qoc <<- function(...) queryObj$cCounts(...)
                inOutFunc$qip <<- function(...) queryObj$populatedRows(...)
                inOutFunc$qop <<- function(...) queryObj$populatedColumns(...)
            }
        }

        if (!is.na(outputDim[2])) {
            ## We can make some callbacks for the mapping matrix
            if (outputDim[2] == 1L) {
                ## Output is represented by rows
                inOutFunc$min <<- function(...) mapObj$cNames(...)
                inOutFunc$mon <<- function(...) mapObj$rNames(...)
                inOutFunc$mic <<- function(...) mapObj$cCounts(...)
                inOutFunc$moc <<- function(...) mapObj$rCounts(...)
                inOutFunc$mip <<- function(...) mapObj$populatedColumns(...)
                inOutFunc$mop <<- function(...) mapObj$populatedRows(...)
            } else {
                ## Output is represented by columns
                inOutFunc$min <<- function(...) mapObj$rNames(...)
                inOutFunc$mon <<- function(...) mapObj$cNames(...)
                inOutFunc$mic <<- function(...) mapObj$rCounts(...)
                inOutFunc$moc <<- function(...) mapObj$cCounts(...)
                inOutFunc$mip <<- function(...) mapObj$populatedRows(...)
                inOutFunc$mop <<- function(...) mapObj$populatedColumns(...)
            }
        }


        if (!is.na(outputDim[3])) {
            ## We can make some callbacks for the ontology matrix
            if (outputDim[3] == 1L) {
                ## Output is represented by rows
                inOutFunc$oin <<- function(...) onotObj$cNames(...)
                inOutFunc$oon <<- function(...) onotObj$rNames(...)
                inOutFunc$oic <<- function(...) onotObj$cCounts(...)
                inOutFunc$ooc <<- function(...) onotObj$rCounts(...)
                inOutFunc$oip <<- function(...) onotObj$populatedColumns(...)
                inOutFunc$oop <<- function(...) onotObj$populatedRows(...)
            } else {
                ## Output is represented by columns
                inOutFunc$oin <<- function(...) onotObj$rNames(...)
                inOutFunc$oon <<- function(...) onotObj$cNames(...)
                inOutFunc$oic <<- function(...) onotObj$rCounts(...)
                inOutFunc$ooc <<- function(...) onotObj$cCounts(...)
                inOutFunc$oip <<- function(...) onotObj$populatedRows(...)
                inOutFunc$oop <<- function(...) onotObj$populatedColumns(...)
            }
        }
        
        outputDim
    },

    weightMatrix = function () {
        "Mutliple-voting weight matrix for mapping matrices, otherwise identity matrix"
        ## This is really just the mapWeights sparse matrix with the
        ## added special sauce of checking to see if the underlying
        ## source matrices have been altered, and if so updating
        ## $mapWeights.
        .updateDerivedStructures()
        mapWeights
    },

    queryToOntology = function () {
        "A weight matrix directly from query IDs to ontology terms"
        ## As with weightMatrix above, this simply returns the
        ## $queryOnto field, but only after updating it if needed.
        .updateDerivedStructures()
        queryOnto
    },

    .noteDiscarded = function (idtype, ids, reason) {
        "Track the reason some IDs are 'discarded' from the analysis"
        alreadyDiscarded <- names(discardedIDs[[ idtype ]])
        newlyDiscarded   <- setdiff(ids, names(alreadyDiscarded))
        nnd <- length(newlyDiscarded)
        if (nnd != 0)  discardedIDs[[ idtype ]] <<-
                           c(discardedIDs[[ idtype ]],
                             stats::setNames(rep(reason, nnd), newlyDiscarded))
    },

    .updateDerivedStructures = function () {
        "Update structures generated from Query, Ontology and optionally Mapping matrix if any of those have changed"
        
        ## Do we have a mapping matrix?
        hasMap <- CatMisc::is.def(mapObj)
        ## Don't do anything if the underlying matrices are unchanged
        ## since the last request.
        if (hasMap) {
            ## We have three matrices: Query, Mapping, Ontology
            if (modState[1] == queryObj$modState &&
                modState[2] == mapObj$modState &&
                modState[3] == ontoObj$modState) return(FALSE)
        } else {
            ## No mapping matrix, just Query and Ontology
            if (modState[1] == queryObj$modState &&
                modState[3] == ontoObj$modState) return(FALSE)
        }
        
### We need to generate the mapping weights. The rows of this matrix
### will be "output" IDs from the Query. The columns will be "input"
### IDs for the Ontology. Both dimensions will be interpreted as
### precisely defining the "worlds" of their relative namespaces.

        qwDef <- CatMisc::is.def(queryWorld)
        if (qwDef) {
            ## If the query world has been defined by the user, we
            ## will ALWAYS accept that as setting the world of IDs.
            activeQueries <<- .standardizeId(queryWorld)
            ## This could result in some "dead end" inputs from the
            ## query (map input rows that connect to zero output
            ## columns). That's ok.
        }
        if (!hasMap) {
            ## No Mapping makes it easy. We just need an identity
            ## matrix representing the world. What world are we using?
            if (!qwDef) {
                ## In the absence of both a user-defined query world
                ## and a mapping matrix, we will presume that the
                ## ontology is defining the query world. Take the
                ## "input" names from the ontology:
                activeQueries <<- .standardizeId( inOutFunc$oin() )
            }
            ## Hand build a simple identity matrix:
            mwLen <- length(activeQueries)
            ijInd <- seq_len(mwLen)
            mapWeights <<- Matrix::sparseMatrix(i=ijInd, j=ijInd, x=1,
                                   giveCsparse=FALSE,
                                   dimnames=list(queryID=activeQueries,
                                                 ontoID=activeQueries))
        } else {
            ## A mapping matrix exists, to convert Query IDs to
            ## Ontology IDs. This is where the main value of
            ## `mapWeights` exists, to numerically manage counting
            ## from one namespace to the other.

            ## What input names are represented in the mapping matrix?
            qids <- .standardizeId( inOutFunc$min() )
            if (!qwDef) {
                ## When a user-defined query world has not been
                ## provided but a mapping matrix has, we will presume
                ## that the mapping matrix is defining the world of
                ## known/utilized query IDs in the "input" dimension:
                activeQueries <<- .standardizeId( qids )
            }
            ## Now begin building a boolean matrix representing the
            ## "valid" in->out / row->col conversions the Mapping
            ## matrix is defining.

            mm <- mapObj$matObj() != 0 # Converts to logical lgTMatrix
            unusedMapRows <- setdiff(qids, activeQueries)
            if (length(unusedMapRows) > 0) {
                ## There are mapping input IDs that are not
                ## represented in our query world. Remove them from
                ## the weight matrix using matched row indices:
                discard <- match(unusedMapRows, qids)
                mm <- mm[ -discard, ,  drop=FALSE ]
                ## Make note of discarded IDs lost from mapping matrix
                .noteDiscarded('queryID', unusedMapRows,
                    "ID in mapping matrix, but not in active queries")
            }

            ## Remove any output IDs (Ontology IDs) that have no
            ## counterparts in the input (colsums are zero).
            populatedColumns <- Matrix::colSums(mm) > 0
            if (!all(populatedColumns)) {
                ## At least one removal
                discard <- colnames(mm)[ !populatedColumns ]
                .noteDiscarded('ontologyID', discard,
                               "Output ID in mapping matrix with no paths from Input IDs")
                mm <- mm[ , populatedColumns, drop=FALSE]
            }
            ## The columns now represent the Ontology world that is
            ## defined/supported by the Query world.

            ## The matrix at this point represents valid paths from
            ## in->out. We now also need to apply weights. These are
            ## simply fractional counts based on the 'multiplicity' of
            ## each row.

            rs <- Matrix::rowSums(mm)
            rowWeights <- ifelse(rs == 0, 0, 1 / rs)
            ## To aid in the matrix math, we'll just make a diagonal
            ## matrix of these weights:
            weightDiag <- Matrix::.sparseDiagonal(length(rowWeights),
                                                  rowWeights)
            ## Replace dimension names so they carry through crossprod:
            sidWD <- .standardizeId(rownames(mm))
            dimnames(weightDiag) <- list(queryID=sidWD,queryID=sidWD)

            mm <- Matrix::crossprod(weightDiag, mm)
            
            ## Finally, we need to add in any activeQueries that are
            ## not represented. They will have rowsums of zero:
            extraQueries <- setdiff(activeQueries, rownames(mm))
            if (length(extraQueries) > 0) {
                nc  <- ncol(mm) # Number of columns in matrix
                nr  <- length(extraQueries) # Number of new rows
                add <- Matrix::sparseMatrix(i=integer(), j=integer(),
                                            x=numeric(),
                                            dims=c(nr, nc),
                                            dimnames=list(queryID=extraQueries,
                                                          ontoID=colnames(mm)))
                ## Stitch the empty rows onto the end of of our matrix
                mm <- Matrix::rbind2(mm, add)
            }
            ## Standardize column names
            colnames(mm) <- .standardizeId(colnames(mm))
            ## Normalize the matrix to dgTMatrix, if it's not already.
            mapWeights <<- as(mm, "dgTMatrix")
        }

        ## mapWeights matrix is normalized so rows are "input" (from
        ## query matrix) and columns are "output" (to ontology
        ## matrix). The row and column names have been normalized with
        ## .standardizeId()

        oids <- colnames(mapWeights) # vector of output/ontology names

        ## Make simplified boolean ontology matrix to use for
        ## calculations. Begin by logical-izing it:
        om <- ontoObj$matrixRaw != 0

        if (outputDim[3] == 1L) {
            ## Transpose the ontology matrix has input in rows, output
            ## in columns:
            om <- Matrix::t( om )
        }
        ## Standardize rownames
        oids2 <- rownames(om) <- .standardizeId(rownames(om))

        ## Only keep rows that are present in the mapping matrix:
        unusedOntoRows <- setdiff(oids2, oids)
        if (length(unusedOntoRows) != 0) {
            ## There are some object IDs defined in the ontology that
            ## we "do not have access to" - they are not represented
            ## in the mapping matrix. Remove them, make a note.
            discard <- match(unusedOntoRows, oids2)
            om <- om[ -discard, ,  drop=FALSE ]
            ## Make note of discarded IDs lost from mapping matrix
            .noteDiscarded('ontologyID', unusedOntoRows,
                           "ID in ontology matrix with no entry in mapping matrix")
        }
        ## Pad out to have the same dimensions as the mapping matrix,
        ## if needed - add in any mapping output IDs that are not in
        ## the simplified ontology matrix:
        extraOntoIds <- setdiff(oids, rownames(om))
        if (length(extraOntoIds > 0)) {
            nc  <- ncol(om) # Number of columns in matrix
            nr  <- length(extraOntoIds) # Number of new rows
            ## Will be lgCMatrix
            add <- Matrix::sparseMatrix(i=integer(), j=integer(),
                                        x=logical(),
                                        dims=c(nr, nc),
                                        dimnames=list(queryID=extraOntoIds,
                                                      ontoID=colnames(om)))
            ## Stitch the empty rows onto the end of of our matrix
            om <- Matrix::rbind2(om, add)
        }
        ## Now fully align/order the ontology matrix rows with the
        ## mapping matrix columns, and store this as the $ontoUse
        ## field. Also standardize to lgTMatrix
        ontoUse <<- as(om[ oids, ,  drop=FALSE ], "lgTMatrix")

### Going to start counting some things. Here are the variables we're
### interested in tallying
        
        ## N = The size of your selected list
        ## i = The count of IDs in your list for an ontology term
        ## n = The total number of IDs held by that ontology term
        ## W = The total number of IDs in the world (aka 'm + n')

        
        ## Let's make a full product of the Mapping Matrix and the
        ## Ontology Matrix. This will allow us to directly convert
        ## query IDs to counts within an ontology term (i)
        queryOnto <<- Matrix::crossprod(Matrix::t(mapWeights), ontoUse)

        ## We can use that matrix to now calculate the total number of
        ## IDs assigned to each ontology term (in the Ontology ID
        ## namespace). We will use generousRound() to integerize these
        ## counts:
        ontoCount <<- generousRound( Matrix::colSums( queryOnto ) )

        ## Are there any ontology terms that end up having *no* IDs
        ## assigned to them?
        noHits <- ontoCount == 0
        if (any(noHits)) {
            ## Ah. At least one ontology is unsupported. There's
            ## really no reason to keep these in the analysis, they
            ## will only serve to slightly penalize p-values when
            ## multiple testing correction takes place. Let's remove
            ## them:
            lostIDs    <- names(ontoCount)[ noHits ]
            ontoCount <<- ontoCount[ !noHits ]
            queryOnto <<- queryOnto[ , !noHits, drop=FALSE ]
            .noteDiscarded('ontologyTerm', lostIDs,
                           "Ontology terms lacking any IDs assigned to them")
        }
        ## Standardize queryOnto
        queryOnto  <<- as(queryOnto, "dgTMatrix")
 
        ## Make some simple structures:
        ontoNames <<- names( ontoCount ) # All surviving ontology terms
        ontoSize  <<- length(ontoNames)  # Number of surviving terms
        
        ## idCount - The fractional count of the output/ontology ID
        ## space, a numeric value between zero and one. These values
        ## represent the maximal count an ontology ID can obtain. The
        ## sum of all of them is worldSize, below, which represents
        ## the maximal count the entire world can obtain.
        idCount <<- pmin( Matrix::colSums(mapWeights), 1 )

        ## worldSize - An integer sum of the maximum number of
        ## ontology IDs representable by this mapping matrix. This is
        ## the "world", or "W", the full count of IDs that can be
        ## selected. Note that some IDs might only contribute a
        ## fractional count - this is intentional! The final value
        ## here will be rounded to an integer.
        worldSize <<- generousRound( sum( idCount ) )
        
        ## Update our modified state to reflect 'now'
        modState <<- c(queryObj$modState,
                       ifelse(hasMap, mapObj$modState, as.integer(NA)),
                       ontoObj$modState )
        TRUE # An update occurred.
        
### TODO : update fields. It's getting messy with additions and removals
### TODO : I feel like I'm leaving a field out ... 
### NEXT : Start tieing this into HGD calculations
### TODO : Strip out un-needed methods
### TODO : TESTS
        
    },

    processAll = function ( force = FALSE, ... ) {
        filter( force = force )
        processListOfLists( force = force )
        adjustResults( force = force )
        topResults( ... )
    },

    threshold = function( x = NA ) {
        if (!is.na(x)) {
            ## New value being provided
            nv <- .normalizeThreshold( x )
            if (!is.na(nv)) logThresh <<- nv[1]
        }
        logThresh
    },
    .normalizeThreshold = function (x) {
        ## Normalize to -log10
        if (x < 1) {
            ## Presume it is being passed as a p-value
            -log(x, base = 10)
        } else {
            ## Presume it is already a -log10 value
            x
        }
    },
    .denormalizeThreshold = function (x) {
        ## Will primarily be used for displaying p-values
        if (x >= 1) {
            ## Presume it is a -log10 value
            x <- 10 ^ -x
        }
        signif(x, 3)
    },
    
    roundLevel = function ( level = NA ) {
        if (!is.na(level)) {
            ## Don't be picky if an integer numeric value is passed:
            if (is.numeric(level) && floor(level) == level)
                level <- as.integer(level)
            roundUp     <<- level
            pseudoRound <<- ifelse(is.integer(level) && level > 0,
                                   (level-1)/level, 0.5)
        }
        roundUp
    },

    ## A pseudocount is added to the value to allow values of >=
    ## 1/roundCount to get rounded up.
    generousRound = function ( x ) {
        rv <- floor(x + pseudoRound)
        ## Set mode in a way that preserves attributes:
        storage.mode(rv) <- "integer"
        rv
    },

    ## Some matrix operations are perfomed with row or column
    ## names. The names may come from different sources and use
    ## different standardization protocols. At the moment I am
    ## concerned about case (1001_at vs 1001_AT or 1001_At), so for
    ## the moment this just normalizes names to lower case

    .standardizeId = function (x) tolower(x),
    
    show = function (...) { cat( .showText(...) ) },

    .resultText = function( color=TRUE ) {
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        if (is.something(resultRaw)) {
            ## Summarize the number of enrichments found that exceeded
            ## the default threshold
            pval <- threshold()
            pvtx <- .denormalizeThreshold(pval)
            rpv  <- resultRaw[ , , 1, drop = FALSE] # Raw result logPvals
            nr   <- sum(abs(rpv) >= pval) # Total number >= threshold
            ## Number of lists that have at least one term >= threshold
            nrl  <- apply(rpv, 1, function(x) sum(abs(x) >= pval))
            nrl  <- sum(nrl > 0)
            ## Number of terms that have at least one term >= threshold
            nrt  <- apply(rpv, 2, function(x) sum(abs(x) >= pval))
            nrt  <- sum(nrt > 0)
            ## Get the dimension name for terms and lists:
            termTxt <- names(dimnames(ontoObj$matrixRaw))[2]
            if (!CatMisc::is.def(termTxt)) termTxt <- "term"
            listTxt <- names(dimnames(queryObj$matrixRaw))[2]
            if (!CatMisc::is.def(listTxt)) listTxt <- "list"
            ## Build a summary line reporting number of positive results
            if (is.something(resultAdj)) {
                apv <- resultAdj[ , , 1, drop = FALSE]
                na   <- sum(abs(apv) >= pval) # Total number >= threshold
                ## Number of lists that have at least one term >= threshold
                nal  <- apply(apv, 1, function(x) sum(abs(x) >= pval))
                nal  <- sum(nal > 0)
                ## Number of terms that have at least one term >= threshold
                nat  <- apply(apv, 2, function(x) sum(abs(x) >= pval))
                nat  <- sum(nat > 0)
                sprintf("p <= %s : %s results = %d %s x %d %s (Raw: %d=%dx%d)",
                        pvtx, doCol(na, color="blue", bgcolor='yellow'),
                        nal, listTxt, nat, termTxt, nr, nrl, nrt)
            } else {
                ## Only raw results available
                sprintf("p <= %s : %s results = %d %s x %d %s %s",
                        pvtx, doCol(nr, color="blue", bgcolor='yellow'),
                        nrl, listTxt, nrt, termTxt,
                        doCol("[Unadjusted]", "red"))
            }
        } else {
            doCol("Not yet analyzed", "purple")
        }
    },

    .showText = function ( compact = FALSE, compactChild = TRUE, color=TRUE ) {
        ## Variable name for use in sample methods
        ## if (!raw) filter()
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        objName <- .self$.selfVarName("myAnalysis")
        objPad  <- "    " # Just for indenting child objects
        name <- param("name")
        msg <- doCol(sprintf("%s analysis results\n",class(.self)[1]),"blue")
        msg <- sprintf("%s    Name: \"%s\"\n", msg, doCol(name, "white"))
        msg <- sprintf("%s    %s\n", msg, .resultText(color=color))
        
        msg <- paste(c(msg, doCol("<#> Query Matrix\n", "blue"),
                           queryObj$matrixText(
                               pad = objPad,
                               compact = compactChild, color=color,
                               fallbackVar = paste(objName,'queryObj',sep='$'))))
        if (CatMisc::is.def(mapObj)) {
            msg <- paste(c(msg, doCol("<#> ID Mapping Matrix\n", "blue"),
                           mapObj$matrixText(pad=objPad, useObj=mapObj$matObj(),
                                            compact=compactChild, color=color,
                                            fallbackVar=paste(objName,'mapObj',sep='$'))))
        }
        msg <- paste(c(msg, doCol("<#> Ontology Matrix\n", "blue"),
                       ontoObj$matrixText(pad = objPad, useObj=ontoObj$matObj(),
                                         compact = compactChild, color=color,
                                         fallbackVar = paste(objName,'ontoObj',sep='$'))))
        filters  <- character()

        if (FALSE) {

### BROKEN. Need to come up with new summary text to report the
### individual Matrix-level filters that are present on the Query,
### Mapping and Ontology filters.
            
        ## Count unique values: https://stackoverflow.com/a/4215196
        doneFilt <- if ("metric" %in% colnames(filterLog)) {
            table(filterLog[ , metric ], filterLog[ , type ])
        } else { character() }
        if (CatMisc::is.def(mapObj)) filters <-
            c(filters,
              .filterHumanText("Minimum mapping score: %s",
                               "minMapMatch", doneFilt, color=color))
        minOntFact <- .ontologyFilterText(color=color)
        if (!is.null(minOntFact)) minOntFact <- paste("  Removing ",
                                                      minOntFact)
        
        filters <-
            c(filters,
              .filterHumanText("Minimum ontology assignment score: %s",
                               "minOntoMatch", doneFilt,color=color),
              minOntFact,
              .filterHumanText("Ontology must not exceed %s%% of world",
                               "maxSetPerc", doneFilt,color=color),
              .filterHumanText("Ontology must have at least %s members",
                               "minSetSize", doneFilt,color=color),
              .filterHumanText("Query must have at least %s ontology terms",
                               "minOntoSize", doneFilt,color=color),
              .filterHumanText2("Recursive filter of mapping matrix",
                                "UnMappedID", doneFilt,color=color),
              .filterHumanText2("IDs with insufficient fractional mapping",
                                "WeakQueryID", doneFilt,color=color)
              )
        }
        
        if (length(filters) != 0) msg <- paste(c(msg, doCol("Filters applied to matrices:\n", "red"), sprintf("  %s\n", filters)))
        worldStats <- character()
        if (is.something(okInput)) worldStats <-
            c(worldStats, paste("Valid Query IDs after filters =",
                                doCol(length(okInput), "red")))
        if (is.something(worldSize)) worldStats <-
            c(worldStats, paste("Total size of world after filters (m+n) =",
                                doCol(worldSize, "red")))
        if (is.something(ontoSize)) worldStats <-
            c(worldStats, paste("Ontology terms surviving filters =",
                                doCol(ontoSize, "red")))
        if (CatMisc::is.def(mapObj) && is.something(roundUp)) {
            worldStats <-
                c(worldStats, paste("Fractional mapped values >=",
                                    doCol(sprintf("1/%d",roundUp), "red"),
                                    "will round up"))
        }
        if (length(worldStats) != 0) msg <- paste(c(msg, doCol("World statistics:\n", "red"), sprintf("  %s\n", worldStats)))

        msg <- paste(c(msg, sprintf("%s: %s$usage()\n",
                       doCol("More options", "red"),
                       doCol(objName, "white"))))
        msg
    },
    .ontologyFilterText = function ( color=TRUE ) {
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        ## human-readable text describing the nature of the ontology
        ## assignment filter
        if (!CatMisc::is.def(param("minOntoMatch"))) return( NULL )
        ## An ontology assignment filter is set
        mom <- param("minOntoMatch")[1]
        lvl <- ontoObj$levels
        if (is.something(lvl)) {
            ## The values of the assignment correspond to levels (and
            ## the presence of a minimum value presumes they are
            ## ordered from "poor" to "good")
            lkept <- seq_len(length(lvl)) >= mom
            fids  <- lvl[ !lkept ]
            if (sum(lkept) / length(lvl) < 0.5) {
                ## Keeping more factor levels than we are discarding
                ftxt <- "not in"
                fids <- lvl[ lkept ]
            } else {
                ## More are discarded
                ftxt <- "in"
                fids <- lvl[ !lkept ]
            }
            sprintf("categories %s %s", ftxt, paste(doCol(
                fids, "red"), collapse = ', '))
        } else {
            ## The minimum filter is some kind of unspecified score,
            ## bigger being better
            sprintf("score < %s", doCol(mom, "red"))
        }
    },
    .filterHumanText = function (fmt = "", var = "", done = NA, key = NA,
                                 color=TRUE) {
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        ## var is either a parameter or field name
        if (hasParam(var)) {
            val <- param(var)
        } else {
            val <- field(var)
        }
        if (!is.something(val)) return( NULL )
        ## Make the basic human-friendly message:
        fmsg <- sprintf(fmt, doCol(val, "red"))
        ## Arg. Vectors will happily tolerate a non-existant name, but
        ## arrays really do not like them (subscript out of
        ## bounds). So check both done and key for presence in the
        ## contingency table
        if (is.something(done) && var %in% rownames(done)) {
            ## A contingency table of counts is provided
            if (!is.character(key)) key <-
                colnames(done[var,done[var,] > 0, drop=FALSE])
            bits = character()
            for (k in key) {
                if (!k %in% colnames(done)) next # This column not present
                dv <- done[ var, k ]
                if (is.something(dv)) {
                    num <- sprintf("%d %s%s", dv, k, ifelse(dv == 1,'','s'))
                    bits <- c(bits, doCol(num, "yellow"))
                }
            }
            if (length(bits) != 0) fmsg <- sprintf("%s %s %s", fmsg,
                          doCol("Removed:", "white"),
                          paste(bits, collapse = ', '))
        }
        fmsg
    },

    .filterHumanText2 = function (fmsg = "", var = "", done = NA, color=TRUE) {
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        ## var is a filter metric name, the contingency table should
        ## be provided and the metric should be present as a row
        if (!is.something(done) || !var %in% rownames(done)) return( NULL )
        ## We will report on all types affected
        key <- colnames(done[var,done[var,] > 0, drop=FALSE])
        bits = character()
        for (k in key) {
            dv <- done[ var, k ]
            if (is.something(dv)) {
                num <- sprintf("%d %s%s", dv, k, ifelse(dv == 1,'','s'))
                bits <- c(bits, doCol(num, "yellow"))
            }
        }
        if (length(bits) != 0) {
            sprintf("%s %s %s", fmsg,
                    doCol("Removed:", "white"), paste(bits, collapse = ', '))
        } else {
            NULL
        }
    },

    .querySummaryText = function ( pad = "", color=TRUE ) {
        .updateDerivedStructures() # Make sure we have up-to-date information
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        name <- doCol(queryObj$param("name", default = "Query"), "yellow")
        dimNames <- names(dimnames(queryObj$matrixRaw))
        txt <- sprintf("%s%s\n", pad, name )
        
    },

    .ontologySummaryText = function ( pad = "", color=TRUE ) {
        .updateDerivedStructures() # Make sure we have up-to-date information
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        name <- doCol(ontoObj$param("name", default = "Ontology"), "yellow")
        dimNames <- names(dimnames(ontoObj$matrixRaw))
        txt <- sprintf("%s%s\n", pad, name )
        if (is.something(ontoSize)) {
            txt <- sprintf("%s%s  %s x %s (filtered from %d x %d)\n", txt, pad,
                           doCol(paste(worldSize, dimNames[1]), "red"),
                           doCol(paste(ontoSize, dimNames[2]), "red"),
                           nrow(ontoObj$matrixRaw), ncol(ontoObj$matrixRaw)
                           )
            txt <- sprintf("%s%s  %s\n", txt, pad, .resultText(color=color))
        } else {
            ## The ontology has been filtered away!
            txt <- sprintf("%s%s  %s\n", txt, pad, doCol(sprintf(
                "Filters have eliminated entire %d x %d ontology!",
                nrow(ontoObj$matrixRaw), ncol(ontoObj$matrixRaw)),
                           color = "red", bgcolor = "yellow"))
        }
        txt
    },

    fisherExact = function(i,N,n,W = NA) {
        ## i  = The number of things in your selection that are 'interesting'
        ##       = phyper 'q' (underenriched) or 'q+1' (over-enriched)
        ## N  = The size of your selection
        ##       = phyper 'k'
        ## n  = The number of interesting things in the world
        ##       = phyper 'm'
        ## W  = The size of the world
        ##       = phyper 'm+n'
        ## So in phyper nomenclature:
        ##  q = i | q = i-1, m = n, n = W - n, k = N
        
        ## Returns p-value as logBase10, with the sign indicating
        ## over/under enrichment. Initially had option to work in
        ## non-logged values, but so much later logic works better in
        ## log space that it was removed (forcing only logged values).

        ## If worldsize was not provided, use the global value
        W[ is.na(W) ] <- worldSize

        ## Allow both N and W to be single values, since they
        ## generally should be constant for a given list /
        ## ontology. We need to recycle them if so
        taskLen <- length(i)
        if (taskLen > 1) {
            if (length(N) == 1) N <- rep(N, taskLen)
            if (length(W) == 1) W <- rep(W, taskLen)
        }

        rv <- ifelse (
            N < 1 | N >= W,
            ## Either nothing (or negative things) were
            ## picked, or the entire world (or more than the
            ## world) was picked. These cases are never
            ## significant (or messed up, but we won't report
            ## errors)
            0,
            ifelse(
                i / N >= n / W,
                ## Our selection has an "interesting ratio" equal to or
                ## greater than the ratio seen in the world
                ## --> Overenriched = return a positive value
                - phyper(i - 1, n, W - n, N,
                         lower.tail = FALSE, log.p = TRUE) / logEtolog10,

                ## Otherwise ...
                ## --> Underenriched = return a negative value
                ## (and compute the other side of the tail)
                phyper(i, n, W - n, N,
                       lower.tail = TRUE, log.p = TRUE) / logEtolog10 ))

        signif( rv, digits = 3 )
    },

    .normalizeList = function( l ) {
        ## Filter out any query IDs that are NOT part of our trimmed,
        ## weighted matrix.
        if (is.null(attr(l, "validInput"))) {
            ## okInput is already standardardized
            validInput <- intersect(.standardizeId(l), okInput)
            vlen       <- length(validInput)
            attr(validInput, "N") <- vlen
            
            ## Note any input IDs that were filtered out:
            attr(l, "rejected")  <- if (vlen == length(l)) {
                ## None removed, all valid
                NULL
            } else {
                ok <- match(okInput, .standardizeId(l))
                ## !is.na(ok) <- the indices we found (ok names)
                ## -ok <- the other indices (non-ok)
                l[ -ok[!is.na(ok) ] ]
            }
            attr(l, "validInput") <- validInput
        }
        l
    },

    processList = function (l, W = NA, format = 'vector',
        no.logging = FALSE) {
        ## Process a single list of IDs against all ontology terms,
        ## returning a vector of logPV (one entry per term)

        validInput <- attr(l, "validInput")
        if (is.null(validInput)) validInput <-
            attr(.normalizeList(l), "validInput")
        
        ## The world size can be over-ridden here
        WS   <- if (is.something(W)) { W } else { worldSize }

        ## vlen is the number of valid IDs provided. If no mapping
        ## matrix is provided, then this is also the same as
        ## "listLen", the size of the selected set. If a mapping
        ## matrix is being used to map to another namespace, then
        ## listLen might be smaller than vlen (if multiple query IDs
        ## map to the same final target ID).
        vlen <- length(validInput)
        if (vlen == 0) {
            ## If there are no valid IDs, return a vector of all zeros
            ## (nothing significant):
            rv <- stats::setNames(rep(0, length(ontoNames)), ontoNames)
            ## Attributes will be needed later:
            attr(rv,"N") <- 0 # List length (count of selected target IDs)
            attr(rv,"W") <- WS      # Size of world (total target IDs)
            attr(rv, 'rejectedIDs') <- attr(l, "rejected")
            attr(rv, 'error') <- "All query IDs were rejected"
            return( rv )
        }
        
        ## Convert the input to a weighted set of genes, and tally up
        ## the number of IDs (potentially fractional if a mapping
        ## matrix is used) assigned to each ontology term
        if (CatMisc::is.def(mapObj)) {
            ## We need to map Query IDs to fractional counts of IDs
            ## used in the ontology.
            wm   <- weightMatrix() # Full weight matrix
            ## Find indices of original weight matrix that match to our
            ## input; Weight matrix names should already be standardized
            vInd <- match(validInput, rownames(wm))
            vim  <- wm[vInd, , drop = FALSE]
            if (length(vim) == 0 || !CatMisc::is.def(nrow(vim))) {
                bogusList <<- l
                rv <- stats::setNames(rep(0, length(ontoNames)), ontoNames)
                ## Attributes will be needed later:
                attr(rv,"N") <- 0 # List length (count of selected target IDs)
                attr(rv,"W") <- WS      # Size of world (total target IDs)
                attr(rv, 'rejectedIDs') <- attr(l, "rejected")
                attr(rv, 'error') <- "Failed to pivot list through weight matrix"
                return( rv )
            }
            ## Find the (potentially fractional) counts of the input
            ## IDs, expressed as a vector
            validInputCounts <- pmin( colSums( vim ), 1 )
            ## Count up the genes and integerize. This is N, the list length
            listLen    <- generousRound(sum(validInputCounts))
            ## Now project the weighted genes onto the ontologies (i):
            myOntoCount <- generousRound(crossprod(
                ontoMatrix( ), validInputCounts )[,1] )
        } else {
            ## No mapping matrix - N is simply the number of queries
            listLen    <- vlen
            ## Ontology counts can be looked up by summing the columns
            ## for matching rows:
            oInd <- match(validInput, .standardizeId(rownames(ontoMatrix())))
            myOntoCount <- colSums( ontoMatrix( )[ oInd, , drop=FALSE])
        }
        ## Calulate the log'ed p-values:
        logPV <- fisherExact( myOntoCount, listLen, ontoCount, WS )
        
        if (!no.logging) {
            if (CatMisc::is.def(mapObj)) {
                msg <- sprintf("List of %d IDs (mapped from %d)",
                               listLen, vlen)
            } else {
                msg <- sprintf("List of %d IDs", vlen)
            }
            msg <- sprintf("%s analyzed against %d ontology terms",msg,
                           ontoSize)
            actionMessage(msg)
        }
        if (grepl('^l', format)) {
            ## List format
            logPV <- list(logPV = logPV, i = myOntoCount, n = ontoCount,
                          N = listLen, W = WS)
        } else if (grepl('^m', format)) {
            ## Matrix format
            logPV <- matrix(
                c(logPV, myOntoCount, ontoCount), ncol = 3, byrow = FALSE,
                dimnames= list(term = ontoNames, metric = c("logPV", "i", "n")))
            ## List length and world size are global values
        } else {
            ## Keep as simple vector
            ## Add names:
            names(logPV) <- ontoNames
        }
        ## Add a few attributes
        attr(logPV,"N") <- listLen # List length (count of selected target IDs)
        attr(logPV,"W") <- WS      # Size of world (total target IDs)
        attr(logPV, 'rejectedIDs') <- attr(l, "rejected")
        logPV
    },

    processListOfLists = function (lol = NULL, lnames = NULL, results = NULL,
        decreasing = FALSE, force = FALSE, ...) {
        usingMethods("processList")
        isDefault <- FALSE
        if (is.null(lol)) {
            ## No list provided, use the object's stored query
            if (is.null(lnames)) {
                ## If other settings are default, return stored raw
                ## result if it is available
                if (CatMisc::is.def(resultRaw) && !force) return( resultRaw )
                isDefault <- TRUE
            }
            lol <- queryObj$matrix( ... )
        }
        ## If an SFMatrix object is provided, get the underlying matrix
        if (inherits(lol, "AnnotatedMatrix")) lol <- lol$matrix( ... )
        if (inherits(lol, c("matrix", "Matrix"))) {
            ## Turn the matrix into a list of character vectors
            ## Build thresholding logic
            min <- param("minQueryScore")
            max <- param("maxQueryScore")
            scoreFilt <- if (CatMisc::is.def(min)) {
                if (CatMisc::is.def(max)) {
                    ## Both min and max filters
                    function(x) x >= min & x <= max 
                    
                } else {
                    ## Min-only filter
                    function(x) x >= min
                }
            } else if (CatMisc::is.def(max)) {
                ## Max-only filter
                function(x) x <= max
            } else {
                ## No filters
                function(x) x != 0

            }
            ## apply was kind of a pain here, since it is difficult to
            ## control if you get a matrix or list back. So use a loop
            ## to normalize. Sorting is irrelevant for standard "bag of
            ## marbles" hypergeometric analysis. However, it is needed
            ## for ranked list analysis (need to implement).
            numL <- ncol(lol)
            extracted <- list()
            for (i in seq_len(numL)) {
                ## 4. Normalize IDs
                ##    3. Extract row names
                ##       2. Sort by value (direction by 'decreasing' param)
                ##          1. Filter values in matrix
                col <- lol[ lol[,i] != 0, i, drop = TRUE ]
                extracted[[ i ]] <- .normalizeList(
                    names(
                        sort(decreasing = decreasing,
                             x = col[ scoreFilt(col) ],
                             )))
            }
            lol <- stats::setNames(extracted, colnames(lol))
        } else if (inherits(lol, "list")) {
            ## This is already what we want
        } else if (is.vector(lol)) {
            ## Simple vector
            lol <- list( myList = lol )
        } else {
            err("processListsOfLists provided with unexpected input",
                class(lol), fatal = TRUE)
        }
        
        rejList <- NULL
        if (is.null(lnames)) {
            ## All lists should be analyzed
            lnames  <- names(lol)
        } else {
            ## Specific request being made
            passed  <- lnames
            lnames  <- intersect(lnames, names(lol))
            ## See if any non-existant requests were made
            if (length(lnames) != length(passed)) {
                ## Where to store rejected list names?
                rejList <- setdiff(passed, lnames)
                lol <- lol[ lnames ]
            }
        }
        numLists <- length(lnames)
        dmsg <- sprintf("Fisher calculation for %d lists x %d terms",
                        numLists, ontoSize)
        name <- param("name")
        if (is.something(name)) dmsg <- paste(dmsg, colorize(name, "white"))
        dateMessage(dmsg)
        ## Establish empty structures to hold results

        results <- array(data = NA, dim = c(numLists, ontoSize, 3),
                         dimnames = list(listName = names(lol),
                             term = ontoNames, metric = c("logPV","i","n")))
        Nvals <- integer(numLists)
        Wvals <- integer(numLists)
        rej   <- list()
        ## Process each list.
        errors <- list()
        for (i in seq_len(numLists)) {
            slice <- processList( lol[[ i ]], no.logging = TRUE,
                                 format = 'matrix')
            ## Add the output to the results array
            results[i, , ] <- slice
            ## Build the relevant attributes
            Nvals[i] <- attr(slice, "N")
            Wvals[i] <- attr(slice, "W")
            rej[[i]] <- attr(slice, 'rejectedIDs')
            etxt <- attr(slice, 'error')
            if (!is.null(etxt)) {
                if (is.null(errors[[ etxt ]])) errors[[ etxt ]] <- 0
                errors[[ etxt ]] <- errors[[ etxt ]] + 1
            }
        }
        for (etxt in colnames(errors)) {
            err(sprintf("[%4d] %s", errors[[ etxt ]], etxt),
                prefix="[ListError]")
        }
        attr(results, "N") <- Nvals
        attr(results, "W") <- Wvals
        attr(results, 'rejectedIDs') <- rej
        attr(results, 'rejectedLists') <- rejList
        dateMessage("Finished", prefix = "  ")
        if (isDefault) resultRaw <<- results
        results
    },

    adjustResults = function ( res = NA, method = 'BY', force = FALSE ) {
        ## Multiple testing adjustment for vectors with +/- sign of
        ## values indicated enriched/underenriched
        isDefault <- FALSE
        if (!CatMisc::is.def(res)) {
            isDefault <- TRUE
            ## If not defined, get the default raw results generated
            ## from processListOfLists(), or return the default
            ## adjusted results if previously defined
            if (CatMisc::is.def(resultAdj) && !force ) return (resultAdj)
            res <- processListOfLists( force = force )
        }
        ## Make a copy of the input:
        rv <- res
        ## The results contain a dimenstion with the logPV, i-values
        ## and n-values. We just want a slice through logPV
        logPV <- rv[ , , 1, drop = FALSE]
        ## Store the sign of each value (to add back later):
        ovrUndr <- ifelse( logPV < 0, -1, 1)

        rv[,,1] <- signif(       # 6. Strip significant digits down to 3
            ovrUndr *            # 5. Replace the over/under sign token
            -log(                # 4. Convert back to log10
                p.adjust(        # 3. Run adjustment
                    10 ^         # 2. Convert back to p-value (^)
                    -abs(logPV), # 1. Un-sign (abs) logPV
                    method = method), base = 10 ), digits = 3)
        if (isDefault) resultAdj <<- rv
        rv
    },

    topResults = function (res=NA, n=100, by="all", force=FALSE,
        add.meta=TRUE, add.raw=TRUE, add.counts=TRUE, caption.name=NA,
        add.urls=FALSE, term.url=NULL, list.url=NULL,
        keep.indices=FALSE, pval=NULL, keep.underenriched=TRUE, ...) {

        ## '...' used to avoid 'unused argument' errors when this
        ## method is ...-called in the same scope as another
        ## ...-called function

        isDefault <- FALSE
        if (!CatMisc::is.def(res)) {
            isDefault <- TRUE
            ## If no results were provided, use default adjusted
            ## results.
            res <- adjustResults( force = force )
        }
        caption <- ""
        fullData <- FALSE
        if (is.something(n)) {
            caption <- sprintf("Top %d", n)
        } else {
            caption <- "All"
            fullData <- TRUE
        }
        caption <- paste(caption, "results")
        
        ## If not provided, use p <= 0.05 as default pval
        pval <- if (is.null(pval)) {
            threshold() } else { .normalizeThreshold( pval ) }
            
        filtBits <- character()
        threshFunc   <- if (is.na(pval)) {
            ## NA = no threshold. Still reject logPV == 0 values
            function(x) x != 0
        } else {
            filtBits <- sprintf("p <= %.3g", 10 ^ -pval)
            if (keep.underenriched) {
                ## Keep both enriched and underenriched
                function(x) { abs(x) >= pval }
            } else {
                ## Keep enriched only
                function(x) { x >= pval }
                filtBits <- c(filtBits, "no under-enrichment")
            }
        }

        getLimit <- function(x) {
            ## Function to find the n'th best result in a set of values:
            tx <- x[ threshFunc(x) ]
            tNum <- length(tx)
            if (tNum == 0) {
                ## No threshold survivors. Use a giant dummy value
                1e99 # This is a log10 value itself, so the max should be ~330
            } else {
                sorted <- sort(abs(x), decreasing = TRUE)
                if (tNum < n || n == 0) { sorted[tNum] } else { sorted[n] }
            }
        }
        ## Function returning indices of all values at least as good as n'th
        getInds  <- function(x, lim) {
            if (keep.underenriched) {
                which( x != 0 & abs(x) >= lim, arr.ind = TRUE)
            } else {
                which( x != 0 & x >= lim, arr.ind = TRUE)
            }
        }
        if (length(filtBits) != 0) caption <-
            paste(caption, sprintf("(%s)", paste(filtBits, collapse = ', ')))
        
        ## The results contain a dimenstion with the logPV, i-values
        ## and n-values. We just want a slice through logPV
        logPV <- res[ , , 1, drop = FALSE]

        rnames   <- rownames(logPV)
        cnames   <- colnames(logPV)
        if (is.null(cnames)) cnames <- character() # from empty results?
        indices  <- NULL
        keycols  <- c("row", "col", 'z')
        if (nrow(logPV) == 0 || ncol(logPV) == 0) {
            ## No results at all. sapply() does unhelpful things, so
            ## just make an empty matrix for indices
            indices <- matrix(nrow=0,ncol=2)
        } else if (!fullData &&
            grepl("(row|list|quer)", by, ignore.case = TRUE)) {
            ## Top n hits on list-by-list basis (lists are in rows)
            by     <- "Row"
            caption  <- paste(caption, "for each query list")
            nR     <- nrow(logPV)
            rInds  <- seq_len(nR)
            ## Track how many top hits were found for each
            ## row. Because we test with >= nth-value, it is possible
            ## for > n hits to be returned. This is a formal field to
            ## suppress "non-local assignment to non-field names"
            ## warnings.
            workSpace <<- integer(nR)
            ## Get all the column indices for each row
            indices <- sapply(rInds, function(ri) {
                                  row <- logPV[ri, ,1]
                                  cInds <- getInds(row, getLimit(row))
                                  workSpace[ri] <<- length(cInds)
                                  cInds
                              })
            cIndices <- unlist(indices)
            numPass  <- length(cIndices)
            if (numPass == 0) {
                ## Nothing passed
                indices <- matrix(integer(), ncol=3, nrow=0)
            } else {
                rIndices <- unlist(sapply(seq_len(nR), function(ri) {
                                              rep(ri, workSpace[ri]) }))
                indices  <- matrix(c(rIndices, cIndices, rep(1, numPass)),
                                   ncol = 3,
                                   dimnames = list(listName = rnames[rIndices],
                                       index = keycols))
            }
        } else if (!fullData &&
                   grepl("(col|term|onto)", by, ignore.case = TRUE)) {
            ## Top n hits on ontology-by-ontology basis (terms are in rows)
            by     <- "Col"
            caption    <- paste(caption, "for each ontology term")
            nC      <- ncol(logPV)
            cInds   <- seq_len(nC)
            ## Track how many top hits were found for each
            ## column. Because we test with >= nth-value, it is
            ## possible for > n hits to be returned. This is a formal
            ## field to suppress "non-local assignment to non-field
            ## names" warnings.
            workSpace <<- integer(nC)
            ## Get all the row indices for each column
            indices <- sapply(cInds, function(ci) {
                                  col <- logPV[ ,ci, 1]
                                  cInds <- getInds(col, getLimit(col))
                                  workSpace[ci] <<- length(cInds)
                                  cInds
                              })
            ## Row indices are a ragged list of passing rows
            ## Col indices are replicated to match up with rows
            rIndices <- unlist(indices)
            numPass  <- length(rIndices)
            if (numPass == 0) {
                ## Nothing passed
                indices <- matrix(integer(), ncol=3, nrow=0)
            } else {
                cIndices <- unlist(sapply(seq_len(nC), function(ci) {
                                              rep(ci, workSpace[ci]) }))
                indices  <- matrix(c(rIndices, cIndices, rep(1, numPass)),
                                   ncol = 3,
                                   dimnames = list(term = cnames[cIndices],
                                       index = keycols))
            }
        } else {
            ## Top n hits across all results at once
            indices <- getInds(logPV, getLimit(logPV))
            if (!fullData) caption <- paste(caption, "across all tests")
        }
        ## indices is a 2-col matrix containing row,col indices
        ## (list,term) of entries within the larger res array that
        ## passed our threshold.
        ri   <- indices[,1]
        ci   <- indices[,2]
        numI <- nrow(indices)
        data <- data.frame(listName = rnames[ri], term = cnames[ci],
                           logPV = logPV[ indices ], stringsAsFactors=FALSE)
        
        if (add.raw) {
            ## Want to include original un-adjusted log-pValues
            data$RawLogPV <- resultRaw[ , , 1, drop = FALSE][ indices ]
        }
        if (add.counts) {
            ## Want to include i,N,n,W counts
            Nvals <- attr(res, "N")
            if (is.null(Nvals)) NVals <- numeric()
            Wvals <- attr(res, "W")
            if (is.null(Wvals)) WVals <- numeric()
            
            data$i <- as.integer(
                res[matrix(c(ri,ci,rep(2,numI)), ncol=3, byrow=FALSE)])
            data$N <- as.integer(Nvals[ri])
            data$n <- as.integer(
                res[matrix(c(ri,ci,rep(3,numI)), ncol=3, byrow=FALSE)])
            data$W <- as.integer(Wvals[ri])
        }
        if (keep.indices) {
            ## Include the array indices for row (list) and col (term)
            data$row <- ri
            data$col <- ci
        }
        tlKey <- c("listName", "term")
        if (nrow(data) != 0) {
            ## We have at least one result; Check if we should add
            ## addtional columns
            if (add.meta) {
                ## Include metadata columns
                for (typ in tlKey) {
                    ## The original matrix for this class of identifier
                    obj  <- if (typ == "term") { ontoObj } else { queryObj }
                    ids  <- unique(data[[typ]]) # the IDs for the data set
                    if (length(ids) == 0) next
                    md   <- obj$metadata( ids ) # All metadata for those ids
                    if (is.null(md)) next
                    mcol <- setdiff(colnames(md), "id") # Non-id metadata cols
                    for (mc in mcol) {
                        if (all(is.na(md[[ mc ]]))) {
                            ## No metadata for this column with these IDs
                            ## Clear it:
                            md[ , (mc) := NULL ]
                        } else {
                            ## rename the column to prefix it with type
                            ## eg "Description" becomes "term Description"
                            setnames(md, mc, paste(typ, mc))
                        }
                    }
                    if (ncol(md) > 1) {
                        ## At least one metadata column is
                        ## available. Splice it into the data table.
                        data <- merge(data, md, all.x = TRUE,
                                      by.x = typ, by.y = "id" )
                    }
                }
            }
            if (add.urls) {
                ## Request to add derivative URL column(s). This is
                ## normally done in htmltable(), but when aggregate
                ## reports are being generated it needs to be done
                ## here.
                urlsAdded <- addUrlColumn(data, term.url, list.url)
                if (!is.null(urlsAdded)) data <- urlsAdded
            }
        }
        ## Order by p-value
        ord   <- order( abs(data$logPV), decreasing = TRUE)
        data  <- data[ ord, ]
        if (is.na(caption.name)) {
            caption <- paste(caption, "for", param("name"))
        } else {
            ## Explicit caption name provided (generally by parent SF object)
            caption <- paste(caption, "for", caption.name)
        }
        attr(data, "caption") <- caption
        attr(data, "by")      <-  by
        attr(data, "pval")    <-  pval
        lastTop <<- data
        data
    },

    .pivotResults = function( x=NA, n=NULL ) {
        ## Uses plyr to pivot result array into a narrow data frame
        ## Dim1 = listName, Dim2 = term
        if (!CatMisc::is.def(x)) x = adjustResults()
        ## https://stackoverflow.com/a/11141750
        require("plyr")
        rv <- adply(x,c(1, 2))
        
        for (metric in c("N","W")) {
            ## Add global N and W values, if available
            vals <- attr(x, metric)
            if (is.something(vals)) {
                rv[[metric]] <- vals[ rv[[1]] ]
            }
        }
        if (is.something(rv$logPV)) {
            ## Sort by pvalue
            rv <- rv[ order(abs(rv$logPV), decreasing = TRUE), ]
        }
        if (is.something(n)) rv <- rv[seq_len(n), ]
        rv
    },

    htmltable = function( res = NA, ... ) {
        if (!CatMisc::is.def(res)) {
            if (CatMisc::is.def(lastTop)) {
                res <- lastTop
            } else {
                err("Can not create HTML table without provided results or running topResults()")
                return(NA)
            }
        }
        setfisher$htmltable( res=res, ... )
    },

    addUrlColumn = function (res, term.url = NULL, list.url = NULL) {
        if (!is.something(res)) return(NULL)
        urls <- list(term     = c(term.url, ontoObj$param('colUrl')),
                     listName = c(list.url, queryObj$param('colUrl')))
        cN   <- colnames(res)
        lcN  <- tolower(cN)
        changes <- 0
        for (col in names(urls)) {
            trgInd <- match(tolower(col), lcN)
            if (is.na(trgInd[1])) next # target column not in output
            urlCol <- paste(col, "URL")
            srcInd <- match(tolower(urlCol), lcN)
            if (is.na(srcInd)) {
                ## URLs were not already explicitly defined. Can we
                ## generate them with a template URL?
                urlFmt <- urls[[col]][ CatMisc::is.def(urls[[col]]) ]
                if (!CatMisc::is.def(urlFmt) || length(urlFmt) == 0) next # No
                ## Apparently yes!
                urlFmt    <- urlFmt[1]
                ## Add a new column containing the URL
                if (grepl('%d', urlFmt)) {
                    ## The URL format has an integer placeholder rather
                    ## than a string. This is used to indicate that an
                    ## integer should be extracted from the presented ID
                    ints <- vapply(res[[trgInd]], function(x) {
                                       as.integer(.parenRE("(\\d+)", x)) }, 0L)
                    res[[urlCol]] <- sprintf(urlFmt, ints)
                } else {
                    ## The format will be presumed to just have a %s holder
                    ## URLencode is not vectorized, so use vapply
                    res[[urlCol]] <- sprintf(urlFmt, vapply(
                        res[[trgInd]], URLencode, "", reserved=TRUE))
                }
                changes <- changes + 1 # Note that we have added a column
            }
        }
        ## If we have modified the DF, return it
        if (changes > 0) return(res)
        NULL # Otherwise use NULL to prevent gratiuitous DF replication
    },

    usage = function() {
        objName <- .self$.selfVarName("myAnalysis")
        fmtM <- sprintf("%s$%%s(%s)\n  %s\n", colorize(objName, "white"),
                       colorize("%s", "purple"), colorize("%s", "cyan"))
        fmtF <- sprintf("%s$%%s\n  %s\n", colorize(objName, "white"),
                        colorize("%s", "cyan"))
        lines <- c(
            sprintf("%s : show summary of object\n", colorize(objName,"white")),
            sprintf(fmtM, "filter", " force=F ", "Apply filters to matrices"),
            sprintf(fmtF, "filterLog", "data.frame of all filtered objects"),
            sprintf(fmtM, "processListOfLists", " ", "Process a list of lists to generate 'raw' HGD values"),
            sprintf(fmtM, "adjustResults", " method='BY' ", "Apply multiple testing adjustment to raw results"),
            sprintf(fmtM, "topResults", " n=100 ", "Report the most significant filtered results"),
            sprintf(fmtM, "processAll", " force=F ", "Run all of the above at once"),
            sprintf(fmtM, "showParameters", " ", "Show configurable parameters for object")
            )
        cat(paste(lines, collapse=""))
    }

    )

## Some utility / testing functions

filteredGoTerms <- function (ana) {
    goFilt <- merge(ana$filterLog[grepl('^GO', ana$filterLog$id)],
                    ana$ontoObj$matrixMD, by = "id")
    goFilt[, c("id", "Description", "filter"), with=FALSE]
}

filteredGenes <- function (ana) {
    goFilt <- merge(ana$filterLog[grepl('^LOC', ana$filterLog$id)],
                    ana$ontoObj$matrixMD, by = "id")
    goFilt[, c("id", "Symbol", "Description", "filter"), with=FALSE]
}


