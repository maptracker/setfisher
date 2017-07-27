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
#' @import CatMisc
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
        setParamList( params=params, ... )
    },

    param = function (key=NA, val=NA, append=FALSE, default=NA,
                      clobber=TRUE, check.class=NULL, is.scalar=NULL,
                      coerce=TRUE) {
        "\\preformatted{
Gets/Sets a parameter. Parameters:
      key - The name of the parameter. If NA will simply return NA. Key
            names are case insensitive.
      val - Optional new value. If not NA, then the parameter value will
            be set, with behavior modified by some of the flags below
   append - Default FALSE, which will cause the parameter to be set to
            'val'. If TRUE, val will be appended to the end of the current
            vector holding the value for 'key'.
  default - Optional value to return if the value of 'key' is 'not defined',
            which corresponds to NA, NULL or zero-length vectors.
  clobber - Default TRUE, which allows an already-set value to be replaced
            with 'val'. Using FALSE is primarily used for managing default
            settings.
 check.class - Default NULL, which will check the parameter definitions
            and use any class found there. If the value is NA or '', then
            there will be no class check. Otherwise, is.class() will be
            tested with the provided class name against the provided
            'val' to see if it matches. If not,  an error will be reported
            and 'key' will not be set. The value 'percent' will be
            interpreted as 'numeric'.
 is.scalar - Default NULL. If TRUE then only val[1] will be taken. If
            FALSE, then all elements of 'val' to be assigned to 'key'.
            If NULL, then if val is only a single element, but matches
            the pattern  '[^]][.+]', then it will be split as a vector
            - see selfSplittingString()
   coerce - Default TRUE, which will attempt to coerce 'val' to
            'check.class' in the event that check.class is not NA
}"
        if (is.na(key)) return( NA )
        ## Always ignore the case of the key:
        key <- tolower(key)
        if (is.null(val)) {
            ## Clear the field
            paramSet[[key]] <<- val
        } else if (any(!is.na(val))) {
            ## Setting a new value
            if (is.null(is.scalar)) {
                val <- selfSplittingString( val )
            } else if (is.scalar) {
                ## Treat the value as a single scalar, and take the
                ## first non-NA/non-NULL value as the one to use. This
                ## allows a preference-ordered list to be provided.
                val <- val[ !is.na(val) ][1]
            }
            if (clobber || !CatMisc::is.def(paramSet[[ key ]])) {
                ## setting clobber to false prevents the value from
                ## being set if one already exists (used for mananging
                ## default settings)
                
                ## Check type of value if a class check was requested
                ## or is stored in the parameter definitions.
                chk <- if (is.null(check.class)) {
                           paramDef[key, "class"] } else { check.class }
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
                if (append && CatMisc::is.def(paramSet[[ key ]])) {
                    ## append flag will extend an existing value
                    paramSet[[ key ]] <<- c(paramSet[[ key ]], val)
                } else {
                    ## Set the value; NULL can be used to clear the field
                    paramSet[[key]] <<- val
                }
            }
        }
        rv <- NA
        if (hasParam(key) && CatMisc::is.def(paramSet[[key]])) {
            rv <- paramSet[[key]]
        } else if (CatMisc::is.def(default)) {
            ## a default value is specified in the event that a
            ## parameter is not set
            rv <- default
        }
        rv
    },

    allParams = function( ) {
        "\\preformatted{
Returns a vector of all parameter keys, either ones that are set or ones
that have a definition set.
}"
        union(rownames(paramDef), names(paramSet))
    },

    hasParam = function(key=NULL) {
        "\\preformatted{
Returns TRUE if the provided key has been set or is in the defaults
      key - Default NULL, should be the key(s) to check
}"
        ## A parameter always exists if there is a definition for
        ## it. This is done to help avoid falling back to looking for
        ## an object field in some functions.
        if (is.null(key)) return( NA )
        is.element( tolower(key), allParams() )
    },

    paramDefinition = function(key, val=NULL) {
        "\\preformatted{
Returns the definition of the parameter, if provided by the code.
      key - Default NULL, should be the key(s) to check
      val - Optional new value to assign
}"
        lkey <- tolower(key)
        if (!is.null(val)) {
            ## Set new value(s)
            paramDef[lkey, "description"] <<- val
            ## Update key column if the requested keys are new
            paramDef[lkey, "key"] <<- ifelse(is.na(paramDef[lkey, "key"]),
                                             key, paramDef[lkey, "key"])
        }
        ## A "naked" NA is storage mode "logical", which behaves
        ## differently than strings or integers when subsetting (it
        ## will return a vector of NAs the length of the number of
        ## rows in the data.frame). To avoid that, I am using the
        ## function below to assure that each value in 'key' will
        ## return one value in the output
        vapply(lkey, function(x) paramDef[x, "description"][1], "")
    },

    paramClass = function(key, val=NULL) {
        "\\preformatted{
Returns the allowed class of the parameter, if provided by the code.
      key - Default NULL, should be the key(s) to check
      val - Optional new value to assign
}"
        lkey <- tolower(key)
        if (!is.null(val)) {
            ## Set new value(s)
            paramDef[lkey, "class"] <<- val
            ## Update key column if the requested keys are new
            paramDef[lkey, "key"] <<- ifelse(is.na(paramDef[lkey, "key"]),
                                             key, paramDef[lkey, "key"])
        }
        ## As above, cumbersome vapply to defend against boolean NA
        vapply(lkey, function(x) paramDef[x, "class"][1], "")
    },

    paramName = function(key, val=NULL) {
        "\\preformatted{
Given a parameter name, returns the name. This is slightly less silly
than it sounds, since names are handled case-insensitively but can
carry a specific case for pretty-printing.
      key - Default NULL, should be the key(s) to check
      val - Optional new value to assign. Must be a case-insensitive
            match to key
}"
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

    setParamList = function (params=NULL, ...) {
        "\\preformatted{
Set one or more parameters provided by a list object. Parameters:
   params - The list holding the parameters, with the names as key names
      ... - dots will be passed to param() for each key/value pair
}"
        if (!is.list(params)) return(NA)
        for (k in names(params)) {
            param(k, params[[k]], ... )
        }
    },

    defineParameters = function(x) {
        "\\preformatted{
Set the parameter definitions (human descriptions) by text block.
   Format of each line is 'keyName [optionalClass] optional description' eg:
        
myFirstKey [integer] Number of widgets to consider
myOtherKey [character] URL for widget lookup
ThatKey Widget asset key, can be text or numeric
}"
        
        if (!CatMisc::is.def(x)) return(NA)
        lines <- unlist(base::strsplit(x, "[\n\r]+"))
        nl <- length(lines)
        keys <- cls <- desc <- rep("", nl)
        nk <- 0
        for (line in lines) {
            keyDat <- CatMisc::parenRegExp("^\\s*(\\S+)(.*)", line)
            key    <- keyDat[1]
            if (CatMisc::is.def(key)) {
                nk <- nk + 1
                keys[nk] <- key
                clsDat <- CatMisc::parenRegExp("^\\s*\\[\\s*(\\S+)\\s*\\](.*)",
                                               keyDat[2])
                if (CatMisc::is.def(clsDat[1])) {
                    cls[nk] <- tolower(clsDat[1])
                    keyDat[2] <- clsDat[2]
                }
                dscDat <-  CatMisc::parenRegExp("^\\s*(.+?)\\s*$", keyDat[2])
                if (CatMisc::is.def(dscDat[1])) desc[nk] <- dscDat[1]
            }
        }
        nks <- seq_len(nk)
        df  <- data.frame(
            key = keys[nks], class = cls[nks], description = desc[nks],
            stringsAsFactors=FALSE)
        rownames(df) <- tolower(keys[nks])
        paramDef <<- df
    },

    showParameters = function () {
        "\\preformatted{
Display available parameters for the object, along with current values
and definitions. Invisibly returns the same text.
}"
        objName <- .selfVarName()
        defFmt  <- paste(c("\n", rep(" ", nchar(objName)[1]), 
                           colorize("# %s", "yellow")), collapse="")
        fmt <- sprintf("%s$param('%s', %s)%%s\n", colorize(objName, "white"),
                       colorize("%s", "magenta"), colorize("%s", "red"))
        lines <- character()
        for (k in allParams() ) {
            v <- param(k)
            com <- attr(v, 'comment')
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
            lines <- c(lines, sprintf(fmt, k , v, desc))
            if (CatMisc::is.something(com[1])) {
                # There is a comment associated with the value
                lines <- c(lines, sprintf(defFmt, strwrap(com[1])))
            }
        }
        txt <- paste(lines, collapse="")
        cat(txt)
        invisible(txt)
    },

    .selfVarName = function( def = "myObj", fallbackVar = "" ) {
        "\\preformatted{
Tries to extract a 'relevant' variable name for displayed help
messages. For example, if show is invoked on object fooThing, this
method is attempting to find the string 'fooThing', so it can show
directly relevant examples like 'fooThing$setWidth()'
      def - Default 'myObj', the string to use if the method fails
            (finding the real object name is not always possible!)
fallbackVar - Another default object name. To be honest, I forget
            why this was needed, but it was. For... reasons.
}"
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
