
logEtolog10 <- log(10) # Natural to base-10 conversion factor

#' SetFisher Analysis
#'
#' Class object associating the components of an enrichment analysis
#' (query lists, ontology, and optional mapping matrix) with the
#' results
#'
#' @field setfisher Pointer back to the parent SetFisher object
#' @field log The SetFisherLogger object holding log (activity)
#'     entries
#' @field queryObj SetFisherMatrix representing the query list(s)
#' @field queryWorld Character verctor describing all objects in the
#'     query namespace. Will default to all unique rows in the mapping
#'     matrix, or all unique rows in the query matrix if no mapping
#'     matrix is present.
#' @field ontoObj SetFisherMatrix representing the mapping of IDs
#'     (from the mapping matrix if present, otherwise from the query
#'     matrix) to ontology terms
#' @field mapObj Optional SetFisherMatrix mapping IDs in the Query
#'     matrix to those used by the Ontology matrix
#' @field mapWeights Weight matrix used to accomodate multiple voting
#'     from the Query namespace to the Ontology namespace
#' @field queryUse Boolean matrix associating queries with
#'     lists. Intended to allow threshold filters to be applied to
#'     ranked input lists, or to allow queryWorld to restrict Query ID
#'     set
#' @field mapUse Boolean matrix indicating allowed mappings from the
#'     Query ID namespace to the Ontology ID namespace. Computed after
#'     applying minMapMatch to the raw Mapping matrix, and eliminating
#'     rows (query input) and columns (ontology output) that lack
#'     entries. Recursive elimination is also currently applied to
#'     remove target IDs that have 'too small' fractional
#'     represenation - this behavior may be removed.
#' @field ontoUse Boolean matrix indicating allowed associations
#'     between IDs (either input IDs, or mapped IDs if a Mapping
#'     matrix is used) and ontology terms. The matrix is filtered
#'     using minOntoMatch (filters original raw matrix based on matrix
#'     score), minOntoSize (minimum number of assigned terms required
#'     to keep an ID), minSetSize (minimum number of IDs required to
#'     keep a term) and maxSetSize (maximum number of IDs allowed for
#'     a term to be kept)
#' @field resultRaw 3D array holding raw p-values from phyper()
#'     calculations, as well as i and n values for each calculation.
#' @field resultAdj p.adjust() processed values from resultRaw
#' @field lastTop data.frame of the most recent values from
#'     topResults()
#' @field isFiltered Boolean flag tracking if filtering has already
#'     been performed on the matrices
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
#' @import data.table
#'
#' @exportClass SetFisherAnalysis
#' 
#' @include SetFisherParamI.R
#' @include SetFisherLoggerI.R
#' @include SetFisherMatrix.R
#' 
#' @export SetFisherAnalysis
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
                    queryObj      = "SetFisherMatrix",
                    queryWorld    = "character",
                    ontoObj       = "SetFisherMatrix",
                    mapObj        = "ANY", # dgTMatrix
                    mapWeights    = "dgCMatrix", # dgTMatrix??

                    queryUse      = "ANY", # dgTMatrix
                    mapUse        = "ANY", # dgTMatrix
                    ontoUse       = "ANY", # dgTMatrix
                    
                    resultRaw     = "array", # [list x ontology x HGD]
                    resultAdj     = "array", # adjusted results
                    lastTop       = "data.frame", # last topResults()
                    
                    isFiltered    = "logical",

                    logThresh     = "numeric", # -log10 value of threshold
                    
                    ## Fractions >= 1/roundUp will be rounded up. The
                    ## default of 2 means that an input query (eg an
                    ## Affy probe set) that maps to two target IDs (eg
                    ## an Entrez gene) can "count" if only 1 out of 2
                    ## are positive.
                    roundUp        = "integer",
                    ## pseudoRound is just (roundUp - 1) / roundUp
                    pseudoRound    = "numeric",
                    log            = "SetFisherLogger",
                    ## qualities of the analysis or world as a whole:
                    ontoSize       = "integer", # Num of ontology terms
                    okInput        = "character", # 'legal' input IDs
                    ontoNames      = "character", # Names of all filtered ontos
                    idCount        = "numeric", # Fractional count for targ IDs
                    ontoCount      = "integer", # Total IDs for each onto
                    worldSize      = "integer", # Total target IDs in world

                    ## Failed list - used for debugging
                    bogusList      = "character",
                    ## scratch variable, to suppress complaints about
                    ## '<<-' needed in an apply() loop
                    workSpace      = "integer"
                    
                    ),
                contains = c("SetFisherLoggerI", "SetFisherParamI")
                )

SetFisherAnalysis$methods(
    
    initialize = function(setfisher = NULL,
        ontology = NULL, query = NULL, idmap = NULL, queryworld = NA,
        minMapMatch    = NA,
        maxSetPerc     = NA,
        minSetSize     = NA,
        minOntoSize    = NA,
        minOntoMatch   = NA,
        round          = 2L,
        param = NA, ... ) {
        usingMethods("param")
        if (!is.def(setfisher)) stop("SetFisherAnalysis entries should be created from a SetFisher object")
        log <<- setfisher$log # Set here so messaging works
        if (!is.def(query)) log$err(
            "SetFisherAnalysis must define 'query' when created", fatal=TRUE)
        if (!is.def(ontology)) log$err(
            "SetFisherAnalysis must define 'ontology' when created", fatal=TRUE)



        
        .self$.setParamDefs("
Name         [character] Optional name for the analysis
minMapMatch  [numeric] Minimum score required to keep a map matrix assignment
maxSetPerc   [percent] Maximum % of world that can be assigned to a term
minSetSize   [integer] Minimum count of set members allowed in a term
minOntoSize  [integer] Minimum number of terms a set member should have
minOntoMatch [numeric] Minimum score requried to keep an ontology matrix assignment
minQueryScore [numeric] Minimum score requried to keep an ID in a query list
maxQueryScore [numeric] Maximum score allowed to keep an ID in a query list
")
        if (is.def(queryworld)) {
            if (length(queryworld) == 0) {
                err("Query world was set to an empty vector - ignoring")
                queryworld <- NA
            } else {
                queryworld <- .standardizeId(queryworld)
            }
        }
        roundLevel( round )
        ## Set some parameters
        mm <- c(minMapMatch)
        if (is.def(idmap)) {
            ## An optional mapping matrix has been provided
            mapObj      <<- idmap
            mm <- c(mm, idmap$param("minMapMatch"))
        }
        param("minMapMatch", mm)
        ## Ontology-level filters
        param("maxSetPerc", c(maxSetPerc, ontology$param("maxSetPerc")))
        param("minSetSize", as.integer(c(
            minSetSize, ontology$param("minSetSize"))))
        param("minOntoSize", as.integer(c(
            minOntoSize, ontology$param("minOntoSize"))))
        param("minOntoMatch", c(minOntoMatch, ontology$param("minOntoMatch")))
        threshold(0.05)
        callSuper(..., setfisher = setfisher, queryWorld = queryworld,
                  queryObj = query, ontoObj = ontology, isFiltered = FALSE )
    },
    
    query       = function (   ) queryObj, # SetFisherMatrix object
    queryMatrix = function (raw = FALSE) { # Matrix
        if (!raw && is.def(queryUse)) return( queryUse )
        queryObj$matrix(raw = raw)
    },
    map       = function (   ) {
        ## SetFisherMatrix object
        if (is.empty.field(mapObj)) { NA } else { mapObj }
    },
    mapMatrix = function (raw = FALSE) { # Matrix
        if (!raw && is.def(mapUse)) return( mapUse )
        mapObj$matrix(raw = raw)
    },
    onto       = function (   ) ontoObj,  # SetFisherMatrix object
    ontoMatrix = function (raw = FALSE) { # Matrix
        if (!raw && is.def(ontoUse)) return( ontoUse )
        ontoObj$matrix(raw = raw)
    },
    weightMatrix = function () {
        ## Is there a better way to check that a field is not initialized?
        if (is.empty.field(mapWeights)) filter()
        mapWeights
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
    
    filter = function ( force = FALSE  ) {
        ## Should only need to filter once
        if (isFiltered && !force) {
            return(filterLog)
        }
        dateMessage(paste("Applying filters to matrices - ",
                          colorize(param("name"), "white")))

        stop("HACKING CODE HERE - moving to AnnotatedMatrix")
        
        ## Set queryUse to a boolean matrix
        queryUse <<- queryObj$matrixRaw != 0
        .shrinkQueryMatrix("absent in user-defined Query World")

        ## Make a boolean matrix (ontoUse) that flags the mappings
        ## between the query (either as provided, or mapped if a
        ## mapping matrix is used) and the ontology. TRUE cells
        ## indicate that the corresponding ontology term is assigned
        ## to the cooresponding query ID.
        minOntoMatch <- param("minOntoMatch")
        if (is.something(minOntoMatch)) {
            ## Filter ontology mappings by a scoring threshold
            ontoUse     <<- ontoObj$matrixRaw >= minOntoMatch
            ## How many ontology assignments did we remove?
            ontoDropped <- nnzero(ontoObj$matrixRaw) - nnzero(ontoUse)
            csR <- colSums(ontoObj$matrixRaw)
            csU <- colSums(ontoUse)
            rsR <- rowSums(ontoObj$matrixRaw)
            rsU <- rowSums(ontoUse)
            if (ontoDropped > 0) {
                ## Filter is causing some ontology *assignments* to be
                ## lost. This may or may not result in ontology terms
                ## or IDs being lost. This row is for the "edges" removed:
                filterDetails(id = "Ontology Matrix", type = "Assignment",
                              filter = sprintf("Ontology assignment score < %s",
                                  minOntoMatch),
                              note = sprintf("%d assignments", ontoDropped))

                ## These are IDs that previously had at least one
                ## assignment to an ontology but now do not. For
                ## reporting purposes, these are the only counts that
                ## will be reported in summaries
                discardedR <- setdiff(names(rsU[rsU == 0]),
                                      names(rsR[rsR == 0]) )
                if (length(discardedR) != 0)
                    filterDetails(id = discardedR, type = "ID",
                                  filter = sprintf("Ontology assignment score < %s",
                                      minOntoMatch), metric = "minOntoMatch" )
                ## These are ontology terms that previously had an
                ## assignment that previously had at least one
                ## assignment to an ID but now do not.
                discardedC <- setdiff(names(csU[csU == 0]),
                                      names(csR[csR == 0]) )
                if (length(discardedC) != 0)
                    filterDetails(id = discardedC, type = "Term",
                                  filter = sprintf("Ontology assignment score < %s",
                                      minOntoMatch), metric = "minOntoMatch"  )
            }
        } else {
            ## Keep all non-zero assignments
            ontoUse <<- ontoObj$matrixRaw != 0
        }
        .shrinkOntologyMatrix(paste("where",.ontologyFilterText()))

        ## Do we have a mapping matrix?
        if (is.def(mapObj)) {
            ## Yes, the namespaces of our query and analysis space differ
            
            ## The mapping matrix plays an important role here in also
            ## defining the size of the (remapped) world
            minMapMatch <- param("minMapMatch")
            if (is.def(minMapMatch)) {
                ## Request to filter the mappings to a minimum score
                mapUse     <<- mapObj$matrixRaw >= minMapMatch
                mapDropped <- nnzero(mapObj$matrixRaw) - nnzero(mapUse)
                csR <- colSums(mapObj$matrixRaw)
                csU <- colSums(mapUse)
                rsR <- rowSums(mapObj$matrixRaw)
                rsU <- rowSums(mapUse)
                if (mapDropped > 0) {
                    ## The filter caused some *mappings* to be
                    ## lost. This may or may not have caused IDs to
                    ## end up being totally excluded. The row below is
                    ## for the "edges" that are lost
                    filterDetails(id = "Mapping Matrix", type = "Mapping",
                                  filter = sprintf("ID mapping score < %s",
                                      minMapMatch),
                                  note = sprintf("%d mappings", mapDropped) )
                    ## These are input IDs that previously had a
                    ## mapping but now do not:
                    discardedR <- setdiff(names(rsU[rsU == 0]),
                                          names(rsR[rsR == 0]) )
                    if (length(discardedR) != 0)
                        filterDetails(id = discardedR, type = "Query ID",
                                      filter = sprintf("ID mapping score < %s",
                                          minMapMatch), metric = "minMapMatch")
                    ## These are output IDs that previously had a
                    ## mapping but now do not. For reporting purposes,
                    ## these are the only counts that will be reported
                    ## in summaries
                    discardedC <- setdiff(names(csU[csU == 0]),
                                          names(csR[csR == 0]) )
                    if (length(discardedC) != 0)
                        filterDetails(id = discardedC, type = "ID",
                                      filter = sprintf("ID mapping score < %s",
                                          minMapMatch), metric = "minMapMatch")

                }

            } else {
                ## No filter, all inputs will contribut to analysis,
                ## provided that they have at least one output
                ## specified in the mapping matrix (value set to
                ## non-zero). We still need to set mapUse to be a
                ## logical sparse matrix
                mapUse <<- mapObj$matrixRaw != 0
            }
            .shrinkMappingMatrix(paste("where score <", minMapMatch))
            
            mids <- .standardizeId(rownames(mapUse))
            if (is.def(queryWorld)) {
                ## User-defined queryWorld - The mapping matrix may
                ## not be completely contained in it
                commonID  <- intersect(mids, queryWorld)
                numDisc   <- length(mids) - length(commonID)
                if (numDisc > 0) {
                    ## Some of the mapping matrix is outside our Query world
                    keepInd   <- match(commonID, mids)
                    discarded <- mids[-keepInd]
                    filterDetails(id = discarded, type = "Query ID", metric = "AlienQuery",
                                  filter = "Mapping Query ID not represented in QueryWorld")
                    actionMessage(sprintf("Removed %d Query IDs from the Mapping Matrix that are not represented in user-supplied world",
                                      numDisc), prefix = "    ")
                    mapUse <<- mapUse[keepInd, , drop = FALSE]
                    .shrinkMappingMatrix("due to non-overlap with query world")
                }
            } else {
                ## If we have not explicitly defined our Query ID
                ## world, then do so using row names on the mapping
                ## matrix
                queryWorld <<- mids
                actionMessage(sprintf("Query World defined as %d Query IDs present in Mapping Matrix",
                                      length(queryWorld)), prefix = "    ")
                .shrinkQueryMatrix("absent in Mapping-defined Query World")
            }
        } else if (is.def(queryWorld)) {
            ## No mapping matrix, but a user-defined query world. Make
            ## sure the ontology matrix is contained within the query
            ## world
            oids      <- .standardizeId(rownames(ontoUse))
            commonID  <- intersect(oids, queryWorld)
            numDisc   <- length(oids) - length(commonID)
            if (numDisc > 0) {
                ## Some of the ontology matrix is outside our Query world
                keepInd   <- match(commonID, oids)
                discarded <- oids[-keepInd]
                filterDetails(id = discarded, type = "ID", metric = "AlienQuery",
                              filter = "Ontology ID not represented in QueryWorld")
                actionMessage(sprintf("Removed %d IDs from the Ontology Matrix that are not represented in user-supplied world",
                                      numDisc), prefix = "    ")
                ontoUse <<- ontoUse[keepInd, , drop = FALSE]
                .shrinkOntologyMatrix("due to non-overlap with query world")
            }
        } else {
            ## No mapping matrix, and no query world either. Define
            ## the query world as all the IDs observed in the
            ## ontology.
            queryWorld <<- .standardizeId(rownames(ontoObj$matrixRaw))
            actionMessage(sprintf("Query World defined as %d IDs present in Ontology Matrix",
                                  length(queryWorld)), prefix = "    ")
            .shrinkQueryMatrix("absent in Ontology-defined Query World")
        }
        

        ## Many of the rows and columns in the matrices could be
        ## reasonably described as "not useful". For example, everyone
        ## should (?) agree that an ontology term with zero assigned
        ## IDs is not useful, same with an ontology that has ALL IDs
        ## assigned to it. Similarly, ontologies with very few or very
        ## many assignments are also not very useful because:

        ## 1. Sparsely-populated ontologies generally fail to reach
        ## statistically significant levels of enrichment. This
        ## threshold is controlled by minSetSize

        ## 2. Heavily-populated ontologies are typically so general as
        ## to not be scientifically informative, like "cell" or
        ## "enzyme". Controlled by maxSetPerc

        ## Excluding ontology terms will provide a (generally small)
        ## statistical "bonus" when multiple testing correction is
        ## applied (fewer tests to account for).

        ## By the same token, there are query IDs that are not useful
        ## for enrichment analysis. These are ones that have low (or
        ## no) level of annotation. IDs that have fewer than
        ## minOntoSize terms assigned to them are removed.

        ## Removing poorly-annotated IDs can have a MAJOR impact on
        ## statistics. For example, un-annotated probe sets in
        ## transcriptional profiling are generally infrequently
        ## expressed, and can even represent spurious (fictional)
        ## genes. These probesets end up clustering at the tail end of
        ## ranked lists, and have the effect of "concentrating" all
        ## other probesets toward the top of the list. The result is
        ## to increase the apparent enrichment of the "real"
        ## probesets. Removing the under-annotated probesets almost
        ## always has the effect of reducing the statistical
        ## significance of ALL ontology terms.

        ## All these filters are applied repeatedly until no further
        ## reduction is seen in either the ontology matrix or the
        ## mapping matrix. A cycle is used because removing ontology
        ## terms can cause a previously passing gene to now fail, and
        ## vise versa.

        cycleTrimming <- 1
        while (cycleTrimming) {
            cycleTrimming <- 0
            cycleTrimming <- cycleTrimming + .pruneOntologyMatrix()
            cycleTrimming <- cycleTrimming + .pruneMappingMatrix()
        }

        ## The rows represent the final number of annotated IDs (genes)
        orn       <- rownames(ontoUse)
        oids      <- .standardizeId(orn)
        if (is.def(mapObj)) {
            ## Trimmed weight matrix rows represent "ok" query
            ## terms. Will intersect user queries with this to remove
            ## rejected IDs from input.

            ## Paranoid that I am messing up somewhere and misaligning
            ## matrices. Ontology rows should match with map columns,
            ## double check here:
            mcn       <- colnames(mapUse)
            mids      <- .standardizeId(mcn)
            if (!identical(oids, mids)) {
                ## GRRR
                err("Mapping and ontology matrices are not aligned. Fixing...",
                    prefix = "[Code Error]")
                commonID <- intersect(oids, mids)
                if (length(commonID) != length(oids)) {
                    ## The ID content is not the same!
                    err("ID count differs between mapping and ontology! Alarmed, fixing...",
                        prefix = "[Code Error]")
                    orInds <- match(commonID, oids)
                    ontoUse <<- ontoUse[orInds, , drop = FALSE]
                }
                mcInds <- match(commonID, mids)
                mapUse <<- mapUse[ , mcInds, drop = FALSE]
                mids   <- commonID
            }
            mrn       <- rownames(mapUse)
            okInput   <<- setNames( .standardizeId(mrn), mrn )
            ## Numeric vector with the fractional mapped count for each gene
            ## that survived trimming
            idCount <<- setNames( pmin( colSums(mapWeights), 1 ), mids)
        } else {
            ## Without a map matrix then the surviving ontology
            ## defines the acceptable IDs
            okInput <<- setNames( oids, orn )
            ## Set idCounts as "1" (pass-through)
            idCount <<- setNames(rep(1, length(okInput)), okInput)
        }
        ## After fractional summing, what is the total number of genes (n+m)?
        worldSize <<- generousRound( sum( idCount ) )
        ## Note all the surviving ontology terms:
        ontoNames <<- colnames(ontoUse)
        ## How many genes are assigned to each ontology term (n)?
        ontoCount <<- generousRound(crossprod(
            as(ontoUse,"dgTMatrix"), idCount )[ , 1] )
        #ontoCount <<- generousRound(crossprod(
        #    as(ontoUse,"dgTMatrix"), idCount )[,1] )
        ## Trimmed ontology matrix reflect final number of ontologies
        ontoSize  <<- length(ontoCount)

        dateMessage("Finished", prefix = "  ")
        isFiltered <<- TRUE
        ## Not really sure what's best to return here. The log, I guess?
        filterLog
    },

    .shrinkOntologyMatrix = function (why='for some un-named reason') {
        ## Shrink the ontology to remove empty rows and columns
        tossRow <- rowSums(ontoUse) == 0
        tossCol <- colSums(ontoUse) == 0
        if (sum(tossRow) + sum(tossCol) == 0) return(NA)
        dimNames <- names(dimnames(ontoObj$matrixRaw))
        actionMessage(sprintf(
            "Ontology removes %d %s and %d %s %s",
            sum(tossRow), dimNames[1], sum(tossCol), dimNames[2],
            why), prefix = "  ")
        ontoUse <<- ontoUse[ !tossRow, !tossCol, drop = FALSE ]
    },
    .shrinkMappingMatrix = function (why='for some un-named reason') {
        ## Shrink the mapping matrix to remove unpopulated rows and columns
        tossRow <- rowSums(mapUse) == 0
        tossCol <- colSums(mapUse) == 0
        if (sum(tossRow) + sum(tossCol) > 0) {
            mapUse <<- mapUse[ !tossRow, !tossCol, drop = FALSE ]
            dimNames <- names(dimnames(mapObj$matrixRaw))
            actionMessage(sprintf(
                "Mapping filter removes %d %s and %d %s %s",
                sum(tossRow), dimNames[1], sum(tossCol), dimNames[2],
                why), prefix = "  ")
        }
    },

    .shrinkQueryMatrix = function (why='for some un-named reason') {
        ## Shrink the query matrix to conform to the Query World
        if (!is.def(queryWorld)) return(NA) # No query world defined
        ## Filter the input queries to assure they are within the
        ## user-defined world of query IDs
        qrn       <- rownames(queryUse)
        qids      <- .standardizeId(qrn)
        commonID  <- intersect(qids, queryWorld)
        numDisc   <- length(qids) - length(commonID)
        if (numDisc == 0) return(NA) # No changes
        ## We have to remove some query IDs from the query matrix
        keepInd   <- match(commonID, qids)
        discarded <- qids[-keepInd]
        filterDetails(id = discarded, type = "Query ID", metric = "AlienQuery",
                      filter = "Query ID not represented in QueryWorld")
        actionMessage(sprintf("Removed %d Query IDs %s", numDisc, why),
                      prefix = "    ")
        if (keepInd == 0) err("All queries have been removed!",
                              prefix = "[Catastrophic Filter]")
        queryUse <<- queryUse[keepInd, , drop = FALSE]
    },

    .pruneOntologyMatrix = function () {
        ## Filter the ontology matrix in a way that recursively
        ## removes rows and columns. Removal of rows may cause some
        ## columns to no longer meet criteria, and vice versa.
        changes   <- 0
        if (length(ontoUse) == 0) return(changes) # Evaporated the entire thing!
        
        message("Filtering assignment matrix", prefix = "  ", color = "blue")

        if (FALSE) {
            ## We will do this in .pruneMappingMatrix()
            if (is.def(mapObj)) {
                ## Keep only IDs that are present as "output" IDs in the
                ## mapping matrix
                nr      <- nrow(ontoUse)
                orn     <- rownames(ontoUse)
                mcn     <- colnames(mapUse)
                okNames <- intersect(.standardizeId(orn), .standardizeId(mcn))
                numDisc <- nr - length(okNames)
                if (numDisc) {
                    ## We need to prune some rows from the ontology
                    keepInd   <- match(okNames, .standardizeId(orn))
                    discarded <- orn[-keepInd]
                    filterDetails(id = discarded, type = "ID",
                                  metric = "MappingUsed",
                                  filter = "ID not represented in mapping matrix")
                    ontoUse <<- ontoUse[keepInd, , drop = FALSE]
                    actionMessage(sprintf("Removed %d unmapped IDs",numDisc),
                                  prefix = "    ")
                    changes <- changes + numDisc
                }
            }
        }

        ## Current size of our working (post-mapping) set 
        totalSetSize <- nrow(ontoUse)
        ## Maximum *number* of set members an ontology should have:
        maxSetPerc   <- param("maxSetPerc")
        maxSetSize   <- generousRound(totalSetSize * maxSetPerc / 100)
        ## Note that as the ontology is pruned, totalSetSize will
        ## decrease, and the maximum allowed number of terms will also
        ## decrease.
        
        minOntoSize  <- param("minOntoSize")
        minSetSize   <- param("minSetSize")
        keepGoing <- TRUE
        while (keepGoing) {
            ## Fast! ~ 2 seconds for large ontologies
            priorChanges <- changes
            if (is.something(minOntoSize)) {
                ## Keep only rows with an adequate number of ontology terms:
                rs      <- rowSums(ontoUse)
                okRow   <- rs >= minOntoSize
                if (any(!okRow)) {
                    ## Some IDs no longer have enough ontology terms
                    orn       <- rownames(ontoUse)
                    discarded <- orn[ ! okRow ]
                    numDisc   <- length(discarded)
                    filterDetails(id = discarded, type = "ID",
                                  filter = sprintf("ID has < %d ontology terms",
                                      minOntoSize), metric = "minOntoSize")
                    ontoUse <<- ontoUse[ okRow, , drop = FALSE]
                    actionMessage(sprintf("Removed %d IDs with < %d ontology terms",
                                        numDisc, minOntoSize),
                                prefix = "    ")
                    changes <- changes + numDisc
                }
            }
            if (is.def(minSetSize)) {
                ## Keep only columns (terms) with a minimum number of
                ## set members
                cs      <- colSums(ontoUse)
                okCol   <- cs >= minSetSize
                if (any(!okCol)) {
                    ## Some ontologies do not have enough assigned IDs
                    ## to be kept
                    ocn       <- colnames(ontoUse)
                    discarded <- ocn[ ! okCol ]
                    numDisc   <- length(discarded)
                    filterDetails(id = discarded, type = "Term",
                                  filter = sprintf("Ontology has < %d IDs",
                                      minSetSize), metric = "minSetSize")
                    ontoUse <<- ontoUse[ , okCol, drop = FALSE ]
                    actionMessage(sprintf("Removed %d ontologies with < %d IDs",
                                        numDisc, minSetSize),
                                prefix = "    ")
                    changes <- changes + numDisc
                }
            }
            if (is.something(maxSetSize)) {
                ## Ontology terms with "too many" loci are generally
                ## not useful and end up needlessly slowing down
                ## analysis (eg 'enzyme' or 'catalytic activity')
                cs      <- colSums(ontoUse)
                okCol   <- cs <= maxSetSize
                numDisc <- length(cs) - sum(okCol)
                if (numDisc) {
                    ## Some ontologies have too many assigned IDs to
                    ## be kept
                    ocn       <- colnames(ontoUse)
                    discarded <- ocn[ ! okCol ]
                    filterDetails(id = discarded, type = "Term",
                                  filter = sprintf("Ontology covers more than %d%% of the world", maxSetPerc ), metric = "maxSetPerc")
                    ontoUse <<- ontoUse[ , okCol, drop = FALSE ]
                    actionMessage(sprintf("Removed %d ontologies with > %d%% world coverage",
                                          numDisc, maxSetPerc),
                                  prefix = "    ")
                    changes <- changes + numDisc
                }
            }
            ## Keep recursing until there are no more changes
            keepGoing <- changes > priorChanges
        }
        if (length(ontoUse) == 0)
            message("Entire ontology matrix has been eliminated by filters",
                    prefix = "    [!!]", bgcolor = "yellow", color = "blue")
        changes
    },
    
    .pruneMappingMatrix = function () {
        changes   <- 0
        if (!is.def(mapObj)) return( changes ) # No map to alter
        if (length(mapUse) == 0 || length(ontoUse) == 0)
            return(changes) # One or both matrices are completely gone
        message("Filtering mapping matrix", prefix = "  ", color = "blue")
        
        ## Make sure that the mapping matrix columns are aligned with
        ## the ontology matrix rows
        mmCol     <- colnames(mapUse)
        mmStnd    <- .standardizeId(mmCol)
        omRow     <- rownames(ontoUse)
        omStnd    <- .standardizeId(omRow)
        commonID  <- intersect(omStnd, mmStnd)
        keepInd   <- match(commonID, mmStnd)
        discarded <- mmCol[-keepInd]
        numDisc   <- length(discarded)
        if (numDisc > 0) {
            ## We have to remove some destination IDs from the mapping matrix
            filterDetails(id = discarded, type = "ID", metric = "UnMappedID",
                          filter = "Map matrix ID not represented in ontology")
            actionMessage(sprintf("Removed %d IDs found in mapping but not ontology",
                                  numDisc),
                          prefix = "    ")
            changes <- changes + length(discarded)
        }
        ## Regardless if any IDs have been removed, assure that the
        ## mapping columns are aligned with the ontology rows:
        mapUse <<- mapUse[ , keepInd, drop = FALSE ]
        
        ## The mapping matrix should have already been converted to a
        ## logical array by $filter(). We now check if there are query
        ## IDs (rows) that no longer have (or never had) *any* target
        ## IDs (columns)
        nonZero   <- rowSums( mapUse ) > 0
        nr        <- nrow(mapUse) # Count of all rows
        numNZ     <- sum(nonZero)           # Count with at least one target
        numDisc   <- nr - numNZ
        if (numDisc) {
            discarded <- rownames(mapUse)[ !nonZero ]
            
            filterDetails(id = discarded, type = "Query ID",
                          metric = "UnMappedID",
                          filter = "Query ID lacking mapped target")
            mapUse <<- mapUse[ nonZero, , drop = FALSE]
            actionMessage(sprintf("%d Query IDs without target IDs removed",
                                  numDisc),
                          prefix = "    ")
            changes <- changes + numDisc
        }
        if (length(mapUse) == 0) {
            actionMessage("Mapping matrix has been completely filtered out!",
                          prefix = "  [!!]", bgcolor='magenta')
            return(0)
        }
        
        ## Now calculate a weight for each query ID based on the
        ## number of target IDs it has
        targCounts  <- rowSums( mapUse )
        targWeights <- 1 / targCounts
        ## Make a sparse diagonal matrix of the weights
        targDiag    <- .sparseDiagonal(length(targWeights), targWeights)
        ## Project these weights into a new numeric matrix
        mapWeights <<- crossprod(targDiag, mapUse)
        ## Normalize the row and column names
        rownames(mapWeights) <<- .standardizeId(rownames(mapUse))
        colnames(mapWeights) <<- .standardizeId(colnames(mapWeights))


#### TO DO - ISSUE TO CONSIDER

        ## Probably should not remove fractional genes - the point was
        ## to allow a large family to "pool their vote" through the
        ## ultimate sum of counts in the ontology, allowing a single
        ## ontological vote to be cast by a family that was hit by
        ## a single probeset
        
#### END ISSUE
        
        ## After calculating fractional counts, we will need to
        ## integer-ize the number of genes, both in the selected set
        ## and in the world. If we use ceiling(), we run the risk of
        ## again over-emphasizing large families (like TSPY or the
        ## protocadherins), which will once more lead to multiple
        ## voting for ontologies linked to those families.

        ## For this reason, we will use round() in counting. This will
        ## immediately eliminate some loci which never reach the level
        ## of at least "half a gene" (or whatever threshold has been
        ## set with roundLevel() ).
        okCol     <- as.logical(generousRound(colSums(mapWeights)))
        numDisc   <- ncol(mapWeights) - sum(okCol)
        if (numDisc) {
            ## There are some targets (eg genes) that will *never*
            ## have enough query IDs (eg probesets) to sum and round up to at
            ## least one gene. Remove them
            ocn       <- colnames(mapWeights)
            discarded <- ocn[ ! okCol ]
            filterDetails(id = discarded, type = "ID",
                          metric = "WeakQueryID",
                          filter ="Target has insufficient total fractional queries")
            mapUse <<- mapUse[ , okCol, drop = FALSE]
            actionMessage(sprintf("Removed %d target IDs with insufficient query representation",
                                  numDisc),
                          prefix = "    ")
            changes <- changes + numDisc
            ## Note that mapWeights is NOT updated here. It is presumed
            ## that will happen on the following iteration.
        }

        ## Finally, tidy up the ontology matrix to account for any IDs
        ## removed from the mapping matrix. In addition to keeping the
        ## three matrices "aligned", this will be needed to allow
        ## recursive trimming in subsequent itterations.
        
        mmStnd    <- .standardizeId(colnames(mapUse))
        if (length(mmStnd) < length(omStnd)) {
            ## We also have ontology IDs that are not represented in
            ## the mapping - remove them.
            keepInd   <- match(mmStnd, omStnd)
            discarded <- omRow[-keepInd]
            numDisc   <- length(discarded)
            filterDetails(id = discarded, type = "ID", metric = "UnMappedID",
                          filter = "Ontology ID not reachable from map matrix")
            actionMessage(sprintf("Removed %d IDs found in ontology but not mapping",
                                  numDisc),
                          prefix = "    ")
            changes <- changes + length(discarded)
            ontoUse <<- ontoUse[ keepInd, , drop = FALSE ]
        }
        if (length(mapUse) == 0)
            message("Entire mapping matrix has been eliminated by filters",
                    prefix = "    [!!]", bgcolor = "yellow", color = "blue")
        changes
    },
    
    roundLevel = function ( round = NA ) {
        if (!is.na(round)) {
            roundUp     <<- round
            pseudoRound <<- ifelse(is.integer(round) && round > 0,
                                   (round-1)/round, 0)
        }
        roundUp
    },

    ## A pseudocount is added to the value to allow values of >=
    ## 1/roundCount to get rounded up.
    generousRound = function ( x ) {
        rv <- round(x + pseudoRound)
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
    
    filterDetails = function ( id = NA, type = NA, filter = NA, note = NA, metric = NA) {
        if (!is.na(filter)) {
            ## Adding new row
            ## debugMessage(sprintf("filterDetails: type = %s, filter = %s, note = %s, metric = %s", paste(type), paste(filter), paste(note), paste(metric)))
            ## debugMessage(sprintf("filterDetails: id = %s, type = %s, filter = %s, note = %s, metric = %s", paste(id), paste(type), paste(filter), paste(note), paste(metric)))
            if (any(is.null(id)) || length(id) == 0) {
                ## Should not happen! This method should only be
                ## called if a non-zero number of IDs are being
                ## filtered out
                err(sprintf(
                    "Filtering noted for metric='%s' on type='%s', but no IDs provided",
                    if (all(is.null(metric))) { "-NULL-" } else { metric },
                    if (all(is.null(type))) { "-NULL-" } else { type }),
                    prefix = "[CodeError]")
                return(NA)
            }
            row <- data.table(id=id, metric = metric, type = type,
                              filter = filter, note = note, key = "id")
            ## Do not add any entries that are already recorded as
            ## filtered (only note the first exclusion)
            newR <- setdiff(row[["id"]], filterLog[["id"]])
            filterLog <<- rbindlist(list(filterLog, row[newR]), fill = TRUE)
            setkeyv(filterLog, "id")
            filterLog
        }
    },

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
            if (!is.def(termTxt)) termTxt <- "term"
            listTxt <- names(dimnames(queryObj$matrixRaw))[2]
            if (!is.def(listTxt)) listTxt <- "list"
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
                           queryObj$.showText(
                               pad = objPad, useObj = queryUse,
                               compact = compactChild, color=color,
                               fallbackVar = paste(objName,'queryObj',sep='$'))))
        if (is.def(mapObj)) {
            msg <- paste(c(msg, doCol("<#> ID Mapping Matrix\n", "blue"),
                           mapObj$.showText(pad = objPad, useObj = mapUse,
                                            compact = compactChild, color=color,
                                            fallbackVar = paste(objName,'mapUse',sep='$'))))
        }
        msg <- paste(c(msg, doCol("<#> Ontology Matrix\n", "blue"),
                       ontoObj$.showText(pad = objPad, useObj = ontoUse,
                                         compact = compactChild, color=color,
                                         fallbackVar = paste(objName,'ontoObj',sep='$'))))
        filters  <- character()
        ## Count unique values: https://stackoverflow.com/a/4215196
        doneFilt <- if ("metric" %in% colnames(filterLog)) {
            table(filterLog[ , metric ], filterLog[ , type ])
        } else { character() }
        if (is.def(mapObj)) filters <-
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
        if (is.def(mapObj) && is.something(roundUp)) {
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
        if (!is.def(param("minOntoMatch"))) return( NULL )
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
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        name <- doCol(queryObj$param("name", default = "Query"), "yellow")
        dimNames <- names(dimnames(queryObj$matrixRaw))
        txt <- sprintf("%s%s\n", pad, name )
        
    },

    .ontologySummaryText = function ( pad = "", color=TRUE ) {
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
        } else if (isFiltered) {
            ## The ontology has been filtered away!
            txt <- sprintf("%s%s  %s\n", txt, pad, doCol(sprintf(
                "Filters have eliminated entire %d x %d ontology!",
                nrow(ontoObj$matrixRaw), ncol(ontoObj$matrixRaw)),
                           color = "red", bgcolor = "yellow"))
        } else {
            ## Has not yet been filtered
            txt <- sprintf("%s%s  %s x %s (unfiltered)\n", txt, pad,
                           doCol(paste(nrow(ontoObj$matrixRaw), dimNames[1]), "red"),
                           doCol(paste(ncol(ontoObj$matrixRaw), dimNames[2]), "red"))
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
            rv <- setNames(rep(0, length(ontoNames)), ontoNames)
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
        if (is.def(mapObj)) {
            ## We need to map Query IDs to fractional counts of IDs
            ## used in the ontology.
            wm   <- weightMatrix() # Full weight matrix
            ## Find indices of original weight matrix that match to our
            ## input; Weight matrix names should already be standardized
            vInd <- match(validInput, rownames(wm))
            vim  <- wm[vInd, , drop = FALSE]
            if (length(vim) == 0 || !is.def(nrow(vim))) {
                bogusList <<- l
                rv <- setNames(rep(0, length(ontoNames)), ontoNames)
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
            if (is.def(mapObj)) {
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
        if (!isFiltered) filter()
        if (is.null(lol)) {
            ## No list provided, use the object's stored query
            if (is.null(lnames)) {
                ## If other settings are default, return stored raw
                ## result if it is available
                if (is.def(resultRaw) && !force) return( resultRaw )
                isDefault <- TRUE
            }
            lol <- queryMatrix( ... )
        }
        ## If an SFMatrix object is provided, get the underlying matrix
        if (inherits(lol, "SetFisherMatrix")) lol <- lol$matrix( ... )
        if (inherits(lol, c("matrix", "Matrix"))) {
            ## Turn the matrix into a list of character vectors
            ## Build thresholding logic
            min <- param("minQueryScore")
            max <- param("maxQueryScore")
            scoreFilt <- if (is.def(min)) {
                if (is.def(max)) {
                    ## Both min and max filters
                    function(x) x >= min & x <= max 
                    
                } else {
                    ## Min-only filter
                    function(x) x >= min
                }
            } else if (is.def(max)) {
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
            lol <- setNames(extracted, colnames(lol))
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
        if (!is.def(res)) {
            isDefault <- TRUE
            ## If not defined, get the default raw results generated
            ## from processListOfLists(), or return the default
            ## adjusted results if previously defined
            if (is.def(resultAdj) && !force ) return (resultAdj)
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
        if (!is.def(res)) {
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
        if (!is.def(x)) x = adjustResults()
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
        if (!is.def(res)) {
            if (is.def(lastTop)) {
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
                urlFmt <- urls[[col]][ is.def(urls[[col]]) ]
                if (!is.def(urlFmt) || length(urlFmt) == 0) next # No
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


