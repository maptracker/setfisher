
typicalEnsSuffix <- '_gene_ensembl'

#' Make Ensembl Matrix
#'
#' Generates an Ensembl AnnotatedMatrix file
#'
#' @param dataset Default \code{NULL}, which will recover all
#'     available datasets from BioMart. Can alternatively be a vector
#'     of one or more datasets.
#' @param \dots Passed on to BioMart functions
#'
#' @details
#'
#' WORK IN PROGRESS
#'
#' @export

makeEnsemblMatrix <- function(dataset=NULL, ...) {

    vers <- biomartVersion(...)
    if (is.na(vers)) return(NA)
    mart <- biomaRt::useMart( names(vers) )

    gVers  <- biomartGenomeVersions(...)
    cNames <- biomartCommonNames(...)

    ## If the dataset has not been specified, use all available datasets
    if (is.null(dataset)) dsets <- names(gVers)

    for (ds in dsets) {
        ## Allow the use of 'hsapiens' instead of 'hsapiens_gene_ensembl':
        needsSfx <- !grepl(paste0(typicalEnsSuffix,'$'), ds)
        if (sum(needsSfx) > 0) ds[needsSfx] <-
                                   paste0(ds[needsSfx], typicalEnsSuffix)
        gv <- gVers[ ds ]
        cn <- cNames[ ds ]
        if (is.na(gv)) {
            message("Unrecognized dataset: ", gv)
            next
        }

        stop("WORK HERE")
        
    }

    
    
}
