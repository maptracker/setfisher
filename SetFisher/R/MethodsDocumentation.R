## ROxygen provides only rudimentary support for documenting Reference
## Class object methods (as of 2017). It will recognize an unassigned
## character string at the beginning of the method function as a
## description, and will parse the parameters. This class is complex
## enough that I wanted more formalized documentation. The blocks
## below define help topics for both RefClass fields and methods:

## * The method name is being put in @name and/or @aliases. This will
##   likely be interpreted by the help system as a simple function in the
##   namespace, but I don't think that will cause problems. It allows
##   \link{}s to be built between topics.
## * I am also specifying the methods under the @method section. This
##   will likely be perceived as an S3/S4 function in the class. I don't
##   _think_ that's an issue, but it might be.
## * It's effectively impossible to use the @usage section, since R_CMD
##   check becomes very unhappy with attempts to formalize usage in
##   'true' object fasion. Instead, usage is manually smuggled into the
##   @details section in a \preformatted{} block.

## ======================= FIELDS =======================

#' Default Query Object
#'
#' Internal SetFisher field holding the Default Query Matrix
#'
#' @name defQryObj
#'
#' @details
#'
#' Either NULL (technically an empty field object) if no default has
#' been defined, or an AnnotatedMatrix object. Set using
#' \link{defaultQuery}.
#'
#' This field is not intended to be directly used. Attempts to
#' manipulate it will likely lead to unpredictable and undesired
#' behavior. Please access it using the \link{defaultQuery}
#' method.
#'
#' @seealso \link{defaultQuery}
#'
#' @examples
#'
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultQuery( s2e )
#' str( sf$defQryObj )
NULL

#' Default Map Object
#'
#' Internal SetFisher field holding the Default Map Matrix
#'
#' @name defMapObj
#'
#' @details
#'
#' Either NULL (technically an empty field object) if no default has
#' been defined, or an AnnotatedMatrix object. Set using
#' \link{defaultMap}.
#'
#' This field is not intended to be directly used. Attempts to
#' manipulate it will likely lead to unpredictable and undesired
#' behavior. Please access it using the \link{defaultMap} method.
#'
#' @seealso \link{defaultMap}
#'
#' @examples
#'
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultMap( s2e )
#' str( sf$defMapObj )
NULL

#' Default Query World
#'
#' Internal SetFisher field holding the default set of query IDs
#'
#' @name defQueryWorld
#'
#' @details
#'
#' Either NULL (technically an empty field object) if no default has
#' been defined, or a character vector of Query IDs. Set with
#' \link{defaultQueryWorld}.
#'
#' This field is not intended to be directly used. Attempts to
#' manipulate it will likely lead to unpredictable and undesired
#' behavior. Please access it using the \link{defaultQueryWorld}
#' method.
#'
#' @seealso \link{defaultQueryWorld}, \link{defaultQuery}
#'
#' @examples
#'
#' sf  <- SetFisher()
#' world <- c("ABC", "XYZ", "ZF12", "HOX4", "OCT4")
#' sf$defaultQueryWorld( world )
#' str( sf$defQueryWorld )
NULL

#' Matrix Cache
#'
#' Internal SetFisher field holding cached matrices
#'
#' @name matrixCache
#'
#' @details
#'
#' A list structure that holds previously requested matrices to allow
#' more rapid recovery and reuse if they are requested multiple times.
#'
#' This field is not intended to be directly used. Attempts to
#' manipulate it will likely lead to unpredictable and undesired
#' behavior. It is built internally through calls to
#' \link{annotatedMatrix}.
#'
#' @seealso \link{annotatedMatrix}, \link{AnnotatedMatrix}
#'
NULL

#' Analysis Cache
#'
#' Internal SetFisher field holding cached analysis objects
#'
#' @name anaCache
#'
#' @details
#'
#' A list structure that holds previously requested SetFisherAnalysis
#' objects to allow more rapid recovery and reuse if they are
#' requested multiple times.
#'
#' This field is not intended to be directly used. Attempts to
#' manipulate it will likely lead to unpredictable and undesired
#' behavior. It is built internally through calls to
#' \link{analysis}.
#'
#' @seealso \link{analysis}, \link{SetFisherAnalysis}
#'
NULL

#' Last Top Results
#'
#' Internal SetFisher field holding the most recent resutls from $topResults
#'
#' @name lastTop
#'
#' @details
#'
#' A data.frame holding the most recent results from
#' \link{topResults}. Designed to allow "efficient" storage of the
#' results (since this package is a ReferenceClass object) for rapid
#' recovery as needed.
#'
#' @seealso \link{topResults}
#'
NULL



## ======================= METHODS ======================

#' Default Query Matrix
#'
#' Set/Get an optional default Query Matrix
#'
#' @name defaultQuery
#' @method defaultQuery SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$defaultQuery( help=TRUE )
#' 
#' queryMat <- myObject$defaultQuery( x=NULL, ... )
#' }
#'
#' If a default query matrix is set, it will be automatically used in
#' cases where one is needed but none is supplied (generally when
#' creating a \link{SetFisherAnalysis})
#'
#' The query is a set of one or more sets of IDs that will be tested
#' for enrichment in one or more ontologies.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     query matrix, or NULL if none is set. Otherwise, will attempt
#'     to recover the AnnotatedMatrix corresponding to \code{x} and
#'     set it as the default.
#' @param \dots Passed to \link{annotatedMatrix}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NULL} if none set, or an \link{AnnotatedMatrix}
#'
#' @seealso \link{defQryObj}, \link{defaultQueryWorld}, \link{defaultMap},
#'     \link{annotatedMatrix}, \link{AnnotatedMatrix}
#'
#' @importFrom CatMisc is.empty.field methodHelp
#'
#' @examples
#' 
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultQuery( s2e )
#' 
NULL

#' Default Map Matrix
#'
#' Set/Get an optional default Map Matrix
#'
#' @name defaultMap
#' @method defaultMap SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$defaultMap( help=TRUE )
#' 
#' mapMat <- myObject$defaultMap( x=NULL, ... )
#' }
#' 
#' If a default mapping matrix is set, it will be automatically used
#' in cases where one is needed but none is supplied (generally when
#' creating a \link{SetFisherAnalysis})
#' 
#' The map is an optional matrix that converts IDs from one
#' "namespace" to another. It is used to connect queries to ontologies
#' when they do not share a namespace. For example, if your queries
#' are sets of Affymetrix probe sets, but your ontology is defined
#' using Ensembl gene accessions, you will need a map that associates
#' Affy probe set IDs with the appropriate Ensembl Gene IDs.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     map matrix, or NULL if none is set. Otherwise, will attempt
#'     to recover the AnnotatedMatrix corresponding to \code{x} and
#'     set it as the default.
#' @param \dots Passed to \link{annotatedMatrix}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NULL} if none set, or an \link{AnnotatedMatrix}
#'
#' @seealso \link{defMapObj}, \link{defaultQuery},
#'     \link{annotatedMatrix}, \link{AnnotatedMatrix}
#'
#' @importFrom CatMisc is.empty.field methodHelp
#'
#' @examples
#' 
#' sf  <- SetFisher()
#' s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )
#' sf$defaultMap( s2e )
#' 
NULL

#' Default Query World
#'
#' Set/Get an optional default Query World
#'
#' @name defaultQueryWorld
#' @method defaultQueryWorld SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$defaultQueryWorld( help=TRUE )
#' 
#' myWorld <- myObject$defaultQueryWorld( x=NULL )
#' }
#' 
#' If a default query world is set, it will be automatically used
#' in cases where one is needed but none is supplied (generally when
#' creating a \link{SetFisherAnalysis})
#'
#' The query world is used to constrict or expand the collection of
#' \strong{all possible} query IDs. For example, you may be running an
#' analysis on mouse genes, but your experimental proceedure
#' explicitly sampled a smaller subset of specific genes. In this
#' case, it is important to constrain the world to that set of genes
#' to properly calculate enrichment statistics.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     query world, or NULL if none is set. Otherwise, will expect a
#'     character vector of IDs that defines the world of query IDs you
#'     wish to use.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NA} if none set, or a character vector of IDs that
#'     defines the world being analyzed.
#'
#' @seealso \link{defQryWorld}, \link{defaultQuery}
#'
#' @importFrom CatMisc is.empty.field methodHelp

#'
#' @examples
#' 
#' sf  <- SetFisher()
#' world <- c("ABC", "XYZ", "ZF12", "HOX4", "OCT4")
#' sf$defaultQueryWorld( world )
NULL

#' Get AnnotatedMatrix
#'
#' Get an AnnotatedMatrix, either from a cache or de novo
#'
#' @name annotatedMatrix
#' @method annotatedMatrix SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$annotatedMatrix( help=TRUE )
#' 
#' myMat <- myObject$annotatedMatrix( x=NA, ... )
#' }
#' 
#' This is a light-weight wrapper around an \code{AnnotatedMatrix()}
#' call, with some additional functionality that an already-loaded
#' matrix will be recovered from an internal cache rather than
#' reloaded.
#' 
#' @param x Default \code{NULL}, which will simply return the default
#'     query world, or NULL if none is set. Otherwise, will expect a
#'     character vector of IDs that defines the world of query IDs you
#'     wish to use.
#' @param \dots Passed to \link{AnnotatedMatrix}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return \code{NULL} if none set, or a character vector of IDs that
#'     defines the world being analyzed.
#'
#' @seealso \link{AnnotatedMatrix}, \link{matrixCache}
#'
#' @import AnnotatedMatrix
#' @importFrom CatMisc methodHelp
#'
#' @examples
#'
#' \dontrun{
#' sf  <- SetFisher()
#' myMat <- sf$annotatedMatix("/data/matrices/Chickens_vs_Cows.mtx")
#' }
#' 
NULL

#' Get SetFisher Analysis
#'
#' Get or create a SetFisher Analysis object
#'
#' @name analysis
#' @method analysis SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$analysis( help=TRUE )
#' 
#' ana <- myObject$analysis(ontology, query=NULL, idmap=NULL,
#'                          name=NULL, queryworld=NULL, ... ) 
#' }
#' 
#' This function is the primary method within SetFisher. It creates a
#' new analysis object, or recovers one that was already created.
#'
#' Note that this function does not actually do any of the analysis,
#' it just mediats the combination of pieces needed to make one (A
#' query matrix, an ontology matrix, and optionally a mapping matrix)
#' and returns an object that will itself be used to analyze for set
#' enrichment.
#' 
#' @param ontology Required, an identifier specifying an Ontology (a
#'     matrix defining one or more gene sets). The identifier can be a
#'     path to the AnnotatedMatrix file, or the name of the ontology
#'     (if it has already been loaded). If ontology is an integer, it
#'     is instead interpreted to be an index value for the internal
#'     analysis cache.
#' @param query Default \code{NULL}, in which case the value defined
#'     by \code{defaultQuery} (if any) will be used. Otherwise and
#'     AnnotatedMatrix (as recovered by \link{annotatedMatrix}) must
#'     be provided.
#' @param idmap Default \code{NULL}, an optional ID mapping
#'     matrix. When NULL, the matrix (if any) stored by
#'     \link{defaultMap} will be used.
#' @param name Default \code{NULL}, an optional human-readable name to
#'     assign to the analysis. If not provided, then a name will be
#'     generated from the name assigned to the ontology.
#' @param queryworld Default \code{NULL}, an optional vector of ID
#'     names to help define (restrict or expand) the "world" that the
#'     query IDs come from. If not provided, will be taken from
#'     \link{defaultQueryWorld}.
#' @param \dots Passed to \link{SetFisherAnalysis}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A \link{SetFisherAnalysis} object.
#'
#' @seealso \link{SetFisherAnalysis}
#'
#' @importFrom CatMisc is.def methodHelp
#'
#' @examples
#'
#' \dontrun{
#' sf  <- SetFisher()
#' myMat <- sf$annotatedMatix("/data/matrices/Chickens_vs_Cows.mtx")
#' }
#' 
NULL

#' Get Analyses from Directory
#'
#' Build analysis objects from all ontology matrices in a directory
#'
#' @name analysesFromDirectory
#' @method analysesFromDirectory SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$analysesFromDirectory( help=TRUE )
#' 
#' anas <- myObject$analysesFromDirectory(dir, query=NULL, idmap=NULL,
#'                                        pattern=NULL, skip=NULL, ...)
#' }
#' 
#' A convienence method that will take pattern-defined matrices from a
#' directory and make an analysis object for each one.
#' 
#' @param dir Required, the path to the directory to process
#' @param query Default NULL, passed on to \link{analysis}
#' @param idmap Default NULL, passed on to \link{analysis}
#' @param pattern Default NULL, which will process all files in
#'     \code{dir}. Otherwise will be passed to \link[base]{list.files}
#' @param skip Default NULL. Can be an optional regular expresion
#'     which when matched against a file name will cause that file to
#'     be skipped.
#' @param \dots Passed to \link{analysis}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A list, with each member being a \code{SetFisherAnalysis}
#'     object.
#'
#' @seealso \link{analysis}
#'
#' @importFrom CatMisc methodHelp
#'
NULL

#' Process All Analyses
#'
#' Filter, process, adjust and report all analyses held by SetFisher
#'
#' @name processAll
#' @method processAll SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$processAll( help=TRUE )
#' 
#' topRes <- myObject$processAll( ... )
#' }
#' 
#' A convienence method that will process all the SetFisherAnalysis
#' objects held by SetFisher, and return a topResults report for them.
#' 
#' @param \dots Passed to \link{topResults}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The output from \link{topResults}
#'
#' @seealso \link{topResults}, \link[SetFisherAnalysis]{filter},
#'     \link[SetFisherAnalysis]{processListOfLists},
#'     \link[SetFisherAnalysis]{adjustResults}
#'
#' @importFrom CatMisc methodHelp
#'
NULL

#' Compact Analysis Names
#'
#' Takes a set of analysis and shortens their names by removing 'common bits'
#'
#' @name analysisNames
#' @method analysisNames SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$analysisNames( help=TRUE )
#' 
#' anaNames <- myObject$analysisNames( trim.names=TRUE, ... )
#' }
#' 
#' While this might appear somewhat frivolous, it's intended to help
#' compact reports that are already rather bulky. If trimming is
#' requested (the default), it's designed to turn name lists like
#' this:
#'
#' \itemize{
#'     \item Human Entrez GeneOntology Analysis
#'     \item Human Entrez MSigDB - Oncogenic Signatures Analysis
#'     \item Human Entrez WikiPathways Analysis
#' }
#'
#' to this:
#' 
#' \itemize{
#'     \item GeneOntology
#'     \item MSigDB - Oncogenic Signatures
#'     \item WikiPathways
#' }
#'
#' @param trim.names Default \code{TRUE}, which will cause the "left"
#'     and "right" sides of all the analysis names to be trimmed to
#'     remove common words. If FALSE, the method isn't too useful.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character vector of analysis names.
#'
#' @seealso \link{dotSharedSubstring}
#'
#' @importFrom CatMisc methodHelp
#'
NULL

#' Top Results Report
#'
#' Generate a data frame summarizing the top results from all analyses
#'
#' @name topResults
#' @method topResults SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$topResults( help=TRUE )
#' 
#' anaNames <- myObject$topResults( ... )
#' }
#' 
#' This method will consider all analysis objects held, and will
#' consolidate their results into a single data.frame. See
#' \link[SetFisherAnalysis]{topResults} (SetFisherAnalysis object) for
#' the recognized parameters that can influence the return value.
#'
#' @param \dots Passed to SetFisherAnalysis'
#'     \link[SetFisherAnalysis]{topResults} method, as well as to
#'     \link{analysisNames}.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A data.frame aggregating the results from all
#'     analyses. This value will also be stored in the \link{lastTop}
#'     field.
#'
#' @seealso \link{lastTop}
#'
#' @importFrom CatMisc methodHelp
#'
NULL

#' TSV Table
#'
#' Create a tab-separated value file from a top results data.frame
#'
#' @name tsvtable
#' @method tsvtable SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$tsvtable( help=TRUE )
#' 
#' tsvPath <- myObject$tsvtable( res=NA, file=NA, force=FALSE, ... )
#' }
#' 
#' Create a simple TSV file from a past analysis.
#'
#' @param res Default \code{NA}, which will use the results in
#'     \link{lastTop}, if available, or call \link{processAll} if
#'     not. Otherwise the user may provide a specific data.frame to
#'     process.
#' @param file Default \code{NA}, which will generate a temporary
#'     file. Otherwise, a specified file path may be provided.
#' @param force Default \code{FALSE}. Passed on to \link{processAll}.
#' @param \dots Passed to \link{processAll} if \code{res} was not
#'     provided.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The path to the TSV file generated.
#'
#' @seealso \link{topResults}
#'
#' @importFrom CatMisc is.something is.def methodHelp
#' @importFrom utils write.table
#'
NULL

#' HTML Table
#'
#' Create a dynamic HTML table from a top results data.frame
#'
#' @name htmltable
#' @method htmltable SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$htmltable( help=TRUE )
#' 
#' htmlPath <- myObject$htmltable(res=NA, file=NA, caption=NULL, show.attr=FALSE,
#'    cols=c("term","term Description", "listName", "logPV", "RawLogPV",
#'           "i","N","i/N","n","W","n/W","Enrich"), ontology.stats=TRUE,
#'           footer=character(), style=character(), force=FALSE,
#'           favicon = "data:Chunk o' Base64 here", ... )
#' }
#' 
#' Create a dynamic HTML table from a data.frame generated by
#' \link{topResults}, decorated with metadata taken from query and
#' ontology matrices.
#'
#' @param res Default \code{NA}, which will use the results in
#'     \link{lastTop}, if available, or call \link{processAll} if
#'     not. Otherwise the user may provide a specific data.frame to
#'     process.
#' @param file Default \code{NA}, which will generate a temporary
#'     file. Otherwise, a specified file path may be provided.
#' @param caption Default \code{NULL}, which will attempt to make a
#'     meaningful caption for the table.
#' @param show.attr Default \code{FALSE}, passed on to
#'     \link[dynamictable]{dynamictable}. The default prevents
#'     data.frame attributes from being printed on the page.
#' @param cols The columns to include in the table.
#' @param ontology.stats Default \code{TRUE}. Mysterious unimplemented
#'     feature. Perhaps I'll remember what I was planning while
#'     driving home.
#' @param footer Default an empty charcter vector. A custom table
#'     footer can be provided if desired. Passed on to
#'     \link[dynamictable]{dynamictable}.
#' @param style Default an empty charcter vector. Custom CSS styles
#'     can be provided if desired. Passed on to
#'     \link[dynamictable]{dynamictable}.
#' @param force Default \code{FALSE}. Passed on to \link{processAll}.
#' @param favicon An explicit \code{data:image/png;base64}
#'     representation of a 16x16 icon to show in the browser tab when
#'     the table is rendered.
#' @param \dots Passed to \link{processAll} if \code{res} was not
#'     provided, and also passed to \link[dynamictable]{dynamictable}.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The path to the HTML file generated.
#'
#' @seealso \link{topResults}
#'
#' @importFrom CatMisc is.something is.def methodHelp
#'
NULL

#' Data Table
#'
#' Create a data.table from a top results data.frame
#'
#' @name datatable
#' @method datatable SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$datatable( help=TRUE )
#' 
#' dt <- myObject$datatable(res=NA, caption=NULL, ... )
#' }
#' 
#' Generate a data.table from a data.frame generated by
#' \link{topResults}, decorated with metadata taken from query and
#' ontology matrices. This uses the
#' \link[SetFisherAnalysis]{datatable} method built into
#' SetFisherAnalysis objects.
#'
#' @param res Default \code{NA}, which will use the results in
#'     \link{lastTop}, if available, or call \link{processAll} if
#'     not. Otherwise the user may provide a specific data.frame to
#'     process.
#' @param caption Default \code{NULL}, which will attempt to make a
#'     meaningful caption for the table.
#' @param \dots Passed to \link{processAll} if \code{res} was not
#'     provided, and also passed to \link[dynamictable]{dynamictable}.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The data.table object
#'
#' @seealso \link{topResults}, \link[SetFisherAnalysis]{datatable}
#'
#' @importFrom CatMisc is.def methodHelp
#'
NULL

#' Fancy Object Text
#'
#' Summarize the object in fancy, human-friendly way
#'
#' @name setFisherText
#' @method setFisherText SetFisher
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$setFisherText( help=TRUE )
#' 
#' dt <- myObject$setFisherText(color=TRUE)
#' }
#' 
#' Create pretty-printed text, optionally colorized, to be used by the
#' RefClass \code{show()} method when summarizing the object.
#'
#' @param color Default \code{TRUE}, which will include ANSI color
#'     tokens in the text for prettier (or messier, depending on your
#'     perspective) text in terminals.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character string.
#'
#' @importFrom CatMisc is.def methodHelp
#'
NULL

## ============ FIELDS - SetFisherAnalysis =============

#' Map Weights Matrix
#'
#' Internal SetFisherAnalysis field, a fractional map from query to ontology IDs
#'
#' @name mapWeights
#'
#' @details
#'
#' This field generally should not be used directly. Instead, using
#' the \link{weightMatrix} method will return this value after
#' checking if it needs to be recalculated.
#' 
#' TODO - Extend documentation!!!
#' 
NULL

#' Query to Ontology Matrix
#'
#' Internal SetFisherAnalysis field, a direct query-to-ontology weighted matrix
#'
#' @name queryOnto
#'
#' @details
#'
#' This field generally should not be used directly. Instead, using
#' the \link{} method will return this value after
#' checking if it needs to be recalculated.
#' 
#' TODO - Extend documentation!!!
#' 
NULL

#' Discarded IDs List
#'
#' Internal SetFisherAnalysis field holding query, mapping and ontology IDs that are discarded / ignored / suppressed
#'
#' @name discardedIDs
#'
#' @details
#'
#' This field notes any identifiers from the Query, Mapping or
#' Ontology matrices that are no longer "available" in the
#' analysis. It is populated by the \link{.noteDiscarded} method,
#' which is called exclusively by \link{.updateDerivedStructures}
#'
#' Identifiers that were previously filtered out using AnnotatedMatrix
#' filter methods will not be included here.
#' 
#' TODO - Extend documentation!!!
#' 
NULL

## ============ METHODS - SetFisherAnalysis =============

#' Generate Weight Matrix
#'
#' Generate the map mapWeights matrix if needed, and then return it
#'
#' @name weightMatrix
#'
#' @details
#'
#' This method will return the \link{mapWeights} matrix if needed, and
#' then return it.
#' 
#' TODO - Extend documentation!!!
#' 
NULL

#' Generate Query To Ontology
#'
#' Generate weighted query-to-ontology matrix if needed, then return it
#'
#' @name queryToOntology
#'
#' @details
#'
#' This method will return the \link{queryOnto} matrix if needed, and
#' then return it.
#' 
#' TODO - Extend documentation!!!
#' 
NULL

#' Record Discarded IDs
#'
#' Internal method that notes IDs discarded during structure processing
#'
#' @name dotNoteDiscarded
#' @aliases .noteDiscarded
#'
#' @details
#'
#' This method will return the \link{queryOnto} matrix if needed, and
#' then return it.
#' 
#' TODO - Extend documentation!!!
#'
#' @keywords internal
NULL

#' Update Derived Structures
#'
#' Internal method that generates derivative matrices and vectors
#'
#' @name dotUpdateDerivedStructures
#' @aliases .updateDerivedStructures
#'
#' @details
#'
#' The two or three matrices used by a SetFisherAnalysis object are
#' pre-processed into a set of derived structures that will be used to
#' perform hypergeometric calculations. This method will create these
#' structures if they do not yet exist, or if changes have been made
#' to the primary matrices (generally by filtering).
#' 
#' TODO - Extend documentation!!!
#'
#' @keywords internal
NULL

