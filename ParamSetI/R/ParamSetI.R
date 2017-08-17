### Generalized paramater storage container

#' Paramter Set Interface
#'
#' Inherited class holding parameter manipulation functions
#'
#' @field paramSet List holding key-value pairs used by an object
#' @field paramDef List holding parameter definitions and class
#'     restrictions
#' @field varName Single string holding extracted variable names. Used
#'     for reporting R code that utilizes actual variable names being
#'     used.
#' 
#' @importFrom CatMisc is.def is.something parenRegExp methodHelp
#' @importFrom methods new setRefClass
#' @importClassesFrom EventLogger EventLogger
#'
#' @examples
#' 
#' ## Demo with inheritance of the class:
#' demo("objectInheritance", package="ParamSetI", ask=FALSE)
#'
#' @export ParamSetI
#' @exportClass ParamSetI
#' 

ParamSetI <-
    setRefClass("ParamSetI",
                fields = list(
                    paramSet = "list",
                    paramDef = "data.frame",
                    varName  = "character" ),
                contains = c("EventLogger") )

ParamSetI$methods(

    initialize = function(params=NA, ... ) {
        "\\preformatted{
Create a new object. Not called directly, but invoked when building an
object with ParamSetI()
     params - Optional list of key/value pairs that will be passed to
              setParamList()
        ... - dots will also be passed on to setParamList()
}"
        callSuper(...)
        paramDef <<- data.frame(key=character(), class=character(),
                                description=character(), stringsAsFactors=FALSE)
        setParamList( params=params, ... )
    },

    param = function (key=NA, val=NA, append=FALSE, default=NA,
                      clobber=TRUE, check.class=NULL, is.scalar=NULL,
                      coerce=TRUE, help=FALSE) {
        "Get or set a parameter"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (is.na(key)) return( NA )
        ## Always ignore the case of the key:
        lkey <- tolower(key)
        if (is.null(val)) {
            ## Clear the field
            paramSet[[lkey]] <<- val
        } else if (any(!is.na(val))) {
            ## Setting a new value
            .addParamKey(key) # Update key column if the requested keys are new
            if (is.null(is.scalar)) {
                val <- selfSplittingString( val )
            } else if (is.scalar) {
                ## Treat the value as a single scalar, and take the
                ## first non-NA/non-NULL value as the one to use. This
                ## allows a preference-ordered list to be provided.
                val <- val[ !is.na(val) ][1]
            }
            if (clobber || !CatMisc::is.def(paramSet[[ lkey ]])) {
                ## setting clobber to false prevents the value from
                ## being set if one already exists (used for mananging
                ## default settings)
                
                ## Check type of value if a class check was requested
                ## or is stored in the parameter definitions.
                chk <- if (is.null(check.class)) {
                           paramDef[lkey, "class"] } else { check.class }
                if (CatMisc::is.something(chk)) {
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
                            suppressWarnings( coerced <-
                                 try({func <- get(sprintf("as.%s", chk))
                                     func(val) }, silent=TRUE) )
                        }
                        if (any(is.na(coerced))) {
                            err(paste("Can not set parameter",key,"to",
                                      val, ": must be of class", chk))
                            return(NA)
                        }
                        ## Success!    
                        val <- coerced
                    }
                }
                if (append && CatMisc::is.def(paramSet[[ lkey ]])) {
                    ## append flag will extend an existing value
                    paramSet[[ lkey ]] <<- c(paramSet[[ lkey ]], val)
                } else {
                    ## Set the value; NULL can be used to clear the field
                    paramSet[[lkey]] <<- val
                }
            }
        }
        rv <- NA
        if (hasParam(lkey) && CatMisc::is.def(paramSet[[lkey]])) {
            rv <- paramSet[[lkey]]
        } else if (CatMisc::is.def(default)) {
            ## a default value is specified in the event that a
            ## parameter is not set
            rv <- default
        }
        rv
    },

    allParams = function( help=FALSE ) {
        "Return all parameter keys, including those that are set or just described"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        union(rownames(paramDef), names(paramSet))
    },

    hasParam = function(key=NULL, help=FALSE) {
        "Test if one or more parameter keys exist in the object"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        ## A parameter always exists if there is a definition for
        ## it. This is done to help avoid falling back to looking for
        ## an object field in some functions.
        if (is.null(key)) return( NA )
        is.element( tolower(key), tolower(allParams()) )
    },

    .addParamKey = function(key) {
        "Make sure that keys are in the definitions table"
        if (!CatMisc::is.something(key)) return(NA)
        lkey <- tolower(key)
        ## Vectorized to allow multiple values in key. If
        ## paramDef$key settings already exist, they will be reset to
        ## themselves. Otherwise the value in argument key will be used
        paramDef[lkey, "key"] <<- ifelse(is.na(paramDef[lkey, "key"]),
                                         key, paramDef[lkey, "key"])
        ## Also update rownames as needed
        isna  <- is.na(rownames(paramDef))
        rownames(paramDef)[isna] <<- tolower(paramDef[isna, "key"])
        sum(isna)
    },

    paramDefinition = function(key, val=NULL, help=FALSE) {
        "Return the definitions, if any, for one or more parameter keys"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        lkey <- tolower(key)
        if (!is.null(val)) {
            ## Set new value(s)
            paramDef[lkey, "description"] <<- val
            .addParamKey(key) # Update key column if the requested keys are new
        }
        ## A "naked" NA is storage mode "logical", which behaves
        ## differently than strings or integers when subsetting (it
        ## will return a vector of NAs the length of the number of
        ## rows in the data.frame). To avoid that, I am using the
        ## function below to assure that each value in 'key' will
        ## return one value in the output
        vapply(lkey, function(x) paramDef[x, "description"][1], "")
    },

    paramClass = function(key, val=NULL, help=FALSE) {
        "Get or set class (storage mode) restrictions of parameters"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        lkey <- tolower(key)
        if (!is.null(val)) {
            ## Set new value(s)
            paramDef[lkey, "class"] <<- val
            .addParamKey(key) # Update key column if the requested keys are new
        }
        ## As above, cumbersome vapply to defend against boolean NA
        vapply(lkey, function(x) paramDef[x, "class"][1], "")
    },

    paramName = function(key, val=NULL, help=FALSE) {
        "Set the capitalization to use for a parameter name"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        lkey <- tolower(key)
        if (!is.null(val)) {
            ## Set new value(s)
            if (tolower(val) == lkey) {
                paramDef[lkey, "key"] <<- val
            } else {
                err(paste("Can not change paramName() of",key,"to",
                          val, ": must only differ in case"))
            }
        }
        ## As above, cumbersome vapply to defend against boolean NA
        vapply(lkey, function(x) paramDef[x, "key"][1], "")
    },

    setParamList = function (params=NULL, help=FALSE, ...) {
        "Set multiple parameters using a list"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!is.list(params)) return(NA)
        for (k in names(params)) {
            param(k, params[[k]], ... )
        }
    },

    defineParameters = function(x, help=FALSE) {
        "Bulk set parameter definitions with a block of text"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!CatMisc::is.def(x)) return(NA)
        lines <- unlist(base::strsplit(x, "[\n\r]+"))
        for (line in lines) {
            keyDat <- CatMisc::parenRegExp("^\\s*(\\S+)(.*)", line)
            key    <- keyDat[1]
            if (CatMisc::is.def(key)) {
                clsDat <- CatMisc::parenRegExp("^\\s*\\[\\s*(\\S+)\\s*\\](.*)",
                                               keyDat[2])
                if (CatMisc::is.def(clsDat[1])) {
                    paramClass(key, clsDat[1])
                    keyDat[2] <- clsDat[2]
                }
                dscDat <-  CatMisc::parenRegExp("^\\s*(.+?)\\s*$", keyDat[2])
                if (CatMisc::is.def(dscDat[1])) paramDefinition(key, dscDat[1])
            }
        }
        paramDef
    },

    showParameters = function ( na.rm=TRUE, help=FALSE ) {
        "Pretty-print parameter names, values and definitions"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        objName <- .selfVarName()
        defFmt  <- paste(c("\n", rep(" ", nchar(objName)[1]), 
                           colorize("# %s", "yellow")), collapse="")
        fmt <- sprintf("%s$param('%s', %s)%%s", colorize(objName, "white"),
                       colorize("%s", "magenta"), colorize("%s", "red"))
        lines <- character()
        for (k in allParams() ) {
            v   <- param(k)
            ## Exclude showing "empty" parameters if requested:
            if (na.rm && (is.null(v) || all(is.na(v)) || all(v == ""))) next
            com <- attr(v, 'comment')[1]
            if (is.character(v)) v <- sprintf("'%s'", v)
            ## Collapse vector values to a single R-like string
            if (length(v) > 1) v <- sprintf("c(%s)", paste(v, collapse=", "))
            ## Has a description been provided for the parameter?
            desc  <- paramDefinition(k)
            desc <- if (CatMisc::is.something(desc)) {
                sprintf(defFmt, desc)
            } else {
                ""
            }
            ## ToDo: Note class restrictions?
            lines <- c(lines, sprintf(fmt, paramName(k) , v, desc))
            if (CatMisc::is.something(com)) {
                # There is a comment associated with the value
                lines <- c(lines, sprintf(defFmt, strwrap(com)))
            }
            lines <- c(lines, "\n")
        }
        txt <- paste(lines, collapse="")
        cat(txt)
        invisible(txt)
    },

    .selfVarName = function( def="myObj", fallbackVar="", help=FALSE ) {
        "Determine the variable name of this object"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!CatMisc::is.something(varName)) {
            ## No attempt to find the object yet. Do so just this once
            for (vn in ls(1)) {
                if (identical(get(vn), .self)) varName <<- vn
            }
            if (!CatMisc::is.something(varName)) {
                ## Still not found
                varName <<- def
            }
        }
        if (CatMisc::is.something(fallbackVar) && varName == def) {
            fallbackVar
        } else {
            varName
        }
    }
)

#' Self-Splitting String
#'
#' Parses a string that defines a split token to turn it into a vector
#'
#' @param x Required, the string vector to proces
#'
#' @details
#'
#' Helper function for ParamSetI, takes a character vector as
#' input. If the input is length 1 and is of format:
#'
#'    [<TOKEN>][<TEXT>]
#'
#' ... then <TEXT> will be split into a vector by <TOKEN>. For example:
#'
#'    [/][x/  y  /z]
#'
#' ... becomes c("x", "  y  ", "z").
#'
#' In all other cases, only a single value will be returned,
#' corresponding to the first non-NA value in x.
#'
#' If the parsed text includes " ## ", then it and the following text
#' will be removed, and the following text will be attached as a
#' 'comment' attribute
#'
#' @return NULL if passed NULL, otherwise a character vector
#'
#' @examples
#'
#' selfSplittingString(" [,][a,b, c ]  ")
#' # "a"   "b"   " c "
#' 
#' ltb <- selfSplittingString("[ and ][lions and tigers and bears] ## oh my")
#' str(ltb, "comment")
#' # "lions" "tigers" "bears"
#' # 'comment' attribute "oh my"
#'
#' # A single return value (vector length 1) is enforced:
#' selfSplittingString(c("a","b","c"))
#' # "a"
#'
#' @importFrom CatMisc parenRegExp is.something
#' 
#' @export

selfSplittingString <- function (x) {
    if (is.null(x)) {
        return( x )
    } else if (length(x) == 1) {
        vecDat <- CatMisc::parenRegExp("^\\s*\\[(.+?)\\]\\[(.+?)\\]\\s*(##\\s+(.+?)\\s*)?$", x)
        if (CatMisc::is.something(vecDat[1])) {
            ## Matches vector pattern
            rv <- unlist(base::strsplit(vecDat[2], vecDat[1], fixed=TRUE))
            if (CatMisc::is.something(vecDat[4])) attr(rv, "comment") <- vecDat[4]
            return( rv )
        }
    }
    rv <- x[ !is.na(x) ][ 1 ]
    hasCom <- CatMisc::parenRegExp("^(.+?)\\s+##\\s+(.+?)\\s*$", rv)
    if (CatMisc::is.something(hasCom[2])) {
        ## Commented 'scalar'
        rv <- hasCom[1]
        attr(rv, "comment") <- hasCom[2]
    }
    rv
}
