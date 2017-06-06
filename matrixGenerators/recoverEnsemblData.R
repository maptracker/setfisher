library("biomaRt")


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


martName <- "ENSEMBL_MART_ENSEMBL"
martList <- listMarts()
martVtxt <- martList[ martList[1] == martName, 2]
if (length(martVtxt) == 0) stop("Could not find Mart with name: ", martName)
martVers <- .parenRE("Ensembl Genes (\\d+)", martVtxt)
if (is.na(martVers)) stop("Expected Mart version of format 'Ensembl Genes ##', got instead '", martVtxt,"'")

martVers <- strtoi(martVers)
geneSfx  <- "_gene_ensembl"

mart <- useMart("ensembl")
dss  <- listDatasets(mart)
numDS    <- nrow(dss)
dss[[1]] <- gsub(geneSfx, '', dss[[1]])    # Remove suffix for compactness
dss[[2]] <- gsub(" genes.+", "", dss[[2]]) # Just keep species part of desc
dss      <- dss[order(toupper(dss[[2]])),] # Sort by species name
rownames(dss) <- seq_len(numDS)
    
.matchDataSet <- function (req) {
    ## Matches a user request for a dataset against those available
    isNum <- .parenRE("^([1-9]\\d*)$", req[1])
    ind   <- NA
    if (!is.na(isNum)) {
        ind <- strtoi(isNum[1])
        if (ind > numDS) {
            message(ind," is larger than the number of available datasets")
            return(NA)
        }
    } else {
        ## A numeric index was not supplied. Check each column for a match
        chkReq <- toupper(req[1])
        chkReq <- gsub(toupper(geneSfx), '', chkReq)
        chkReq <- gsub(toupper(" genes.+"), "", chkReq)
        for (j in seq_len(ncol(dss))) {
            ind <- which(toupper(dss[[j]]) == chkReq)
            if (length(ind) > 1) {
                stop("Multiple matches found for ",req)
            } else if (length(ind) == 1) {
                break
            }
        }
    }
    if (length(ind) == 1 && !is.na(ind)) {
        as.character(dss[ind,])
    } else {
        NA
    }
}

setDat <- NA
if (exists("dataset") && !is.na(dataset)) {
    ## dataset has already been defined in current search space
    setDat <- .matchDataSet(dataset)
    if (is.na(setDat[1])) {
        print(dss)
        stop("Requested dataset '",dataset,"' was not recognized. Please pick from list above")
    }
} else {
    ## Show list of available datasets and let user pick
    .pickDataSet <- function() {
        print(dss)
        readline("Please choose a dataset from the list above
   (any of the four columns can be provided as a selection): ")
    }
    while (is.na(setDat[1])) {
        setDat <- .matchDataSet( .pickDataSet() )
    }
}

message("

WORKING HERE
* Build base file name
* Recover basic gene information
* WHERE ARE THE UNOFFICIAL SYMBOLS ?!?

")
str(setDat)

