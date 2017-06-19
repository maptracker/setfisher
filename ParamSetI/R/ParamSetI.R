### Generalized paramater storage container

#' SetFisher Paramter Interface
#'
#' Inherited class holding parameter manipulation functions
#'
#' @field paramSet List holding key-value pairs used by an object
#' @field paramDef List holding default key-value pairs
#' @field varName Single string holding extracted variable names. Used
#'     for reporting R code that utilizes actual variable names being
#'     used.
#' 
#' @importFrom dynamictable dynamictable
#' @importFrom methods setRefClass new
#'
#' @import data.table
#'

SetFisherParamI <-
    setRefClass("SetFisherParamI",
                fields = list(
                    paramSet = "list",
                    paramDef = "data.frame",
                    varName  = "character"
                    ),
                contains = c("SetFisherLoggerI")

                )


SetFisherParamI$methods(

    param = function (key = NA, val = NA, append = FALSE, default = NA,
                      clobber = TRUE, check.class = NA, is.scalar = TRUE,
                      coerce = TRUE) {
        if (is.na(key)) return( NA )
        ## Always ignore the case of the key:
        key <- tolower(key)
        if (is.null(val)) {
            ## Clear the field
            paramSet[[key]] <<- val
        } else if (any(!is.na(val))) {
            ## Setting a new value
            if (is.scalar) {
                ## Treat the value as a single scalar, and take the
                ## first non-NA/non-NULL value as the one to use. This
                ## allows a preference-ordered list to be provided.
                val <- val[ !is.na(val) ][1]
            }
            if (clobber || !is.def(paramSet[[ key ]])) {
                ## setting clobber to false prevents the value from
                ## being set if one already exists (used for mananging
                ## default settings)
                vecDat <- .parenRE("^\\s*\\[(.+?)\\]\\[(.+?)\\]\\s*$", val)
                if (is.something(vecDat[1])) {
                    ## The construct "[,][A,B,C]" is used to specify a
                    ## vector of values. The first [] block specifies
                    ## the text delimiter, the second the values
                    val <- unlist(base::strsplit(vecDat[2], vecDat[1]))
                }
                ## Check type of value if a class check was requested
                ## or is stored in the parameter definitions.
                chk <- ifelse(is.na(check.class), paramDef[key, "class"],
                              check.class)
                if (is.def(chk)) {
                    ## 'percent' is a conveinence class to help with
                    ## automated value reporting.
                    if (chk == "percent") chk <- "numeric"
                    if (!inherits(val, chk)) {
                        ## Eh. was not what we expected.
                        coerced <- NA
                        if (coerce) {
                            ## User is willing to try forcing it. This
                            ## is the default since many settings will
                            ## have been grepped out of text files and
                            ## will begin life as characters.
                            coerced <- try({func <- get(sprintf("as.%s", chk))
                                func(val) }, silent = TRUE)
                        }
                        if (is.na(coerced)) {
                            err(paste("Can not set parameter",key,"to",
                                      val, ": must be of class", chk))
                            return(NA)
                        }
                        ## Success!    
                        val <- coerced
                    }
                }
                if (append && is.def(paramSet[[ key ]])) {
                    ## append flag will extend an existing value
                    paramSet[[ key ]] <<- c(paramSet[[ key ]], val)
                } else {
                    ## Set the value; NULL can be used to clear the field
                    paramSet[[key]] <<- val
                }
            }
        }
        rv <- NA
        if (hasParam(key) && is.def(paramSet[[key]])) {
            rv <- paramSet[[key]]
        } else if (is.def(default)) {
            ## a default value is specified in the event that a
            ## parameter is not set
            rv <- default
        }
        rv
    },

    hasParam = function(key) {
        ## A parameter always exists if there is a definition for
        ## it. This is done to help avoid falling back to looking for
        ## an object field in some functions.
        tlk <- tolower(key)
        tlk %in% rownames(paramDef) || tlk %in% names(paramSet)
    },

    .setParamList = function (x, ...) {
        ## Set one or more parameters provided by a list object
        if (!is.list(x)) return(NA)
        for (k in names(x)) {
            param(k, x[[k]], ... )
        }
    },

    .setParamDefs = function(x) {
        ## Set the parameter definitions by text block. Format of each
        ## line is "keyName [optionalClass] optional description" eg:
        
        ## myFirstKey [integer] Number of widgets to consider
        ## myOtherKey [character] URL for widget lookup
        ## ThatKey Flexible reporting value for this 'n' that
        
        if (!is.def(x)) return(NA)
        lines <- unlist(base::strsplit(x, "[\n\r]+"))
        nl <- length(lines)
        keys <- cls <- desc <- rep("", nl)
        nk <- 0
        for (line in lines) {
            keyDat <- .parenRE("^\\s*(\\S+)(.*)", line)
            key    <- keyDat[1]
            if (is.def(key)) {
                nk <- nk + 1
                keys[nk] <- key
                clsDat <- .parenRE("^\\s*\\[\\s*(\\S+)\\s*\\](.*)", keyDat[2])
                if (is.def(clsDat[1])) {
                    cls[nk] <- tolower(clsDat[1])
                    keyDat[2] <- clsDat[2]
                }
                dscDat <-  .parenRE("^\\s*(.+?)\\s*$", keyDat[2])
                if (is.def(dscDat[1])) desc[nk] <- dscDat[1]
            }
        }
        nks <- seq_len(nk)
        df  <- data.frame(
            key = keys[nks], class = cls[nks], description = desc[nks],
            stringsAsFactors=FALSE)
        rownames(df) <- tolower(keys[nks])
        paramDef <<- df
    },

    .selfVarName = function( def = "myObj", fallbackVar = "" ) {
        ## In illustrative output I want to include the 'actual'
        ## variable names for the objects being used. That is, instead
        ## of showing "myMatrix$metadata('GeneSymbol')" I want the
        ## code to find the actual variable name the user picked, eg
        ## "hgu133a$metadata('GeneSymbol')". This will make it easier
        ## to explore the code with copy-paste

        ## This turns out to be a REAL pain to do. I spent a lot of
        ## time looking at deparse / sys.call to figure out what the
        ## "real" variable name was. In particular this approach fails
        ## inside implicit calls to show(). If show is explicitly
        ## called - "show(myObject)" - then the variable name can be
        ## found on the call stack. For some reason if myObject was
        ## printed (show()n) implicitly by calling just "myObject",
        ## then the stack "does not go back far enough"

        ## So instead I am going to brute force look at the variable
        ## space. I am only looking in the global environment (1), and
        ## there are many ways this might not be quite right, but it
        ## is close enough. I hope.
        if (!is.something(varName)) {
            ## No attempt to find the object yet. Do so just this once
            for (vn in ls(1)) {
                if (identical(get(vn), .self)) varName <<- vn
            }
            if (!is.something(varName)) {
                ## Still not found
                varName <<- def
            }
        }
        if (is.something(fallbackVar) && varName == def) {
            fallbackVar
        } else {
            varName
        }
    }
)
