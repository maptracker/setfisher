#' Event Logger
#'
#' Utility ReferenceClass object for messaging and recording of
#' events. Use \code{$message()} for general logging (stores event and
#' reports to console), see Methods section for variants.
#'
#' A $show() method has been set, so simply evaluating an EventLogger
#' object on the command line will pretty-print the result. The
#' data.table holding the log information can be directly accessed
#' in field $log.
#'
#' @field log The data.table holding log messages
#' @field useCol Logical flag to indicate if color should be used in
#'     messaging
#' @field vb Logical flag indicating if verbose messaging should be
#'     active
#' @field colMap An internally generated list that maps color names to
#'     crayon functions
#'
#' @importFrom methods new setRefClass
#' @importFrom CatMisc is.def methodHelp
#' @importFrom data.table data.table rbindlist
#' @import crayon
#'
#' @examples
#'
#' myEL <- EventLogger()
#' myEL$message("Did something important")
#' Sys.sleep(3)
#' myEL$actionMessage("Something emphatic has happened")
#' Sys.sleep(1)
#' myEL$dateMessage("Here's a date stamp")
#' myEL$debugMessage("Remember to comment this out in production")
#'
#' # Pretty print the log, including an elapsed time:
#' myEL
#' # Expose the underlying data.table:
#' myEL$log
#'
#' ## Demo with inheritance of the class:
#' demo("objectInheritance", package="EventLogger", ask=FALSE)
#'
#' 
#' @export EventLogger
#' @exportClass EventLogger
#'

EventLogger <-
    setRefClass("EventLogger",
                fields = list(
                    useCol   = "logical",
                    log      = "data.table",
                    vb       = "logical",
                    colMap   = "list",
                    EvLogObj = "ANY",
                    varName  = "character"
                    ),
                )

EventLogger$methods(
    
    initialize = function(useColor=NULL, verbose=NULL, log=NULL, ...) {
        "\\preformatted{
Create a new object using EventLogger():
        log - Optional EventLogger object. This is used if EventLogger is an
              inherited ('contains') class in another RefClass object, and you
              wish that object to share the event log from a previously created
              object.
   useColor - Defaults to TRUE, toggles if crayon::-based coloring is applied
              to printed messages. If 'log' is provided, changes will affect
              that shared object
    verbose - Defaults to TRUE, toggles if log events are also printed to the
              terminal. If 'log' is provided, changes will affect that shared
              object
}"
        EvLogObj <<- log
        if (is.null(EvLogObj)) {
            ## Make a log table from scratch - this will be utilized
            log <<- data.table::data.table(Date=Sys.time(),
              Message="Log initialized", key="Date")
        } else {
            ## Make a "dummy" log table to fill the field; It will NOT
            ## be used, instead the $log in EvLogObj will be utilized
            log <<- data.table::data.table(Date=Sys.time(),
              Message="STUB TABLE - see $EvLogObj for shared table", key="Date")
        }
        ## Manage useColor and verbose flags through methods. Default
        ## is TRUE for both, but can be user set.
        useCol   <<- TRUE
        vb       <<- TRUE
        .self$useColor(useColor)
        .self$verbose(verbose)
        ## Set the string-to-colorizerMethod lookup list:
        colNameToFunc( )
    },

    message = function(msg="No message provided!", prefix=NULL,
                       color=NULL, bgcolor=NULL, datestamp=FALSE,
                       fatal=FALSE, collapse=" ", help=FALSE) {
        
        "Display a formatted message, and store it in the log table"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!is.null(prefix)) prefix <- c(prefix, " ") # Optional prefix
        ti <- Sys.time()
        ds <- if (datestamp) {
            ## Append date stamp to message
            ds <- colorize(as.character(ti), "yellow")
            c(" - ", ds)
        } else {
            NULL
        }
        m <- colorize(c(prefix, paste(msg, collapse=collapse)), color)
        if (datestamp) m <- c(m, colorize(ds, "yellow"))
        m <- colorize(m, NULL, bgcolor)
        if (fatal) {
            stop(m)
        } else if (verbose()) {
            base::message(m)
        }
        ## Add entry to log table
        if (!is.null(EvLogObj)) {
            ## Use a shared object
            EvLogObj$log <<-
                data.table::rbindlist(list(EvLogObj$log, data.table::data.table(
                 Date=ti, paste(msg, collapse=" "))))
            invisible(EvLogObj$log)
        } else {
            log <<-
                data.table::rbindlist(list(log, data.table::data.table(
                 Date=ti, paste(msg, collapse=" "))))
            invisible(log)
        }
    },

    dateMessage = function ( msg="No message provided!", help=FALSE, ... ) {
        "Calls message() with datestamp=TRUE"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        message(msg=msg, datestamp=TRUE, ...)
    },

    actionMessage = function (msg="No message provided!!", prefix='[+]',
        color="red", help=FALSE, ...) {
        "Calls message with a '[+]' prefix and red coloring"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        message(msg=msg, prefix=prefix, color=color, ...)
    },

    debugMessage = function (msg="No message provided!!", prefix='[DEBUG]',
        color="white", bgcolor="blue", help=FALSE, ...) {
        "Calls message with a '[DEBUG]' prefix and white/blue coloring"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        message(msg=msg, prefix=prefix, color=color, bgcolor=bgcolor, ...)
    },

    err = function (msg="No message provided!!", prefix='[ERROR]', help=FALSE,
    ...) {
        "Calls message with an '[ERROR]' prefix and red/yellow coloring"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        message(msg=msg, prefix=prefix, collapse="\n",
                color="red", bgcolor="yellow", ...)
    },

    colorMap = function(color, bg=FALSE, help=FALSE) {
        "Convert a color name into a crayon colorizing function"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!useColor() || ! CatMisc::is.def(color)) return(NA)
        if (is.function(color[1])) return( color[1] )
        key   <- "FG"
        color <- tolower(color[1])
        if (bg || grepl('^bg', color)) {
            ## Use the background methods
            key   <- "BG"
            ## Lookup keys are standardized to discard the leading bg:
            color <-  gsub('^bg','', color)
        }
        colNameToFunc()[[ key ]][[ color ]]
    },

    colorize = function(msg="", color=NULL, bgcolor=NULL, help=FALSE) {
        "Use crayon to add ANSI color codes to text"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!is.character(msg)) msg <- as.character(msg)
        fgFn <- colorMap(color, FALSE)
        if (is.function(fgFn)) msg <- fgFn(msg)
        bgFn <- colorMap(bgcolor, TRUE)
        if (is.function(bgFn)) msg <- bgFn(msg)
        msg
    },
    
    useColor = function(newval=NULL, help=FALSE) {
        "Get or set the flag determining if messages are colorized"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!is.null(EvLogObj)) {
            ## Use a shared object
            return( EvLogObj$useColor(newval=newval) )
        } else if (!is.null(newval)) {
            nv <- as.logical(newval)[1]
            if (is.na(nv)) {
                err("useColor() should be provided with a boolean argument")
            } else {
                useCol <<- nv
            }
        }
        invisible(useCol)
    },

    verbose = function(newval=NULL, help=FALSE) {
        "Get or set the flag determining if messages are displayed"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        if (!is.null(EvLogObj)) {
            ## Use a shared object
            return( EvLogObj$verbose(newval=newval) )
        } else if (!is.null(newval)) {
            nv <- as.logical(newval)[1]
            if (is.na(nv)) {
                err("verbose() should be provided with a boolean argument")
            } else {
                vb <<- nv
            }
        }
        invisible(vb)
    },

    tidyTime = function (x=NULL, pad=0, help=FALSE) {
        "Reports a time interval with unit management and colorization"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        unit  <- "s"
        color <- "yellow";
        if (is.null(x) || is.na(x)) {
            unit  <- '?'
            x     <- -999
            color <- 'magenta'
        } else if (x < 1) {
            unit  <- "ms"
            x     <- x * 1000
            color <- "green"
        } else if (x > 60) {
            unit  <- "min"
            x     <- x / 60
            color <- "red"
        }
        ## Colorizing includes non-printing tokens, such as '\033[31m',
        ## that still "count" when computing padding. For this reason, pad
        ## the string before we colorize.
        colorize(sprintf(paste('%', pad, 's', sep=""),
                         sprintf("%.3f %s", x, unit)), color)
    },

    showLog = function (help=FALSE, ...) {
        "Pretty-prints the log, including total elapsed time"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        cat( .self$logText(...) )
    },

    show = function (...) {
        "A wrapper for showLog"
        showLog( ... )
    },

    logText = function (width=0.7 * getOption("width"),
                        relative=TRUE, pad=11, n=0, help=FALSE) {
        "Formats the log as a human readable two-column table"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        usingMethods("tidyTime") # Needed for use in apply
        head <- colorize("Activity log:", "blue")
        ## Use a shared log object, if available
        use <- if(is.null(EvLogObj)) { log } else { EvLogObj$log }
        ## Nicely format the log
        HMS <- use$Date
        nl  <- length(HMS)
        tot <- difftime(HMS[nl], HMS[1], units='secs')
        if (relative) {
            if (nl < 2) {
                HMS <- colorize("Empty Log", "magenta")
            } else {
                ## We want to report the time used for each step
                HMS <- difftime(HMS[ 2:nl ], HMS[ 1:(nl-1) ], units="secs")
                ## The times are deltas, so one line will not have a
                ## time associated with it. This seems to "line up"
                ## best if the non-stamped event is the last one,
                ## presuming that messages are "emitted" at the start
                ## of an event, rather than the end.
                HMS <- c(unlist(lapply(HMS, tidyTime, pad = pad)),
                         strrep(" ",pad) )
            }
        } else {
            HMS <- colorize(format(HMS, "%H:%M:%S"), "yellow")
        }
        msgs <- use$Message
        if (n > 0) {
            msgs <- tail(msgs, n)
            HMS  <- tail(HMS, n)
        }
        ## Word wrap and pad the messages
        wrapCol <- paste("\n", strrep(" ", pad + 3), sep="")
        msg <- unlist(lapply(apply(as.matrix(msgs),
                                   1, strwrap, width=width),
                             paste, collapse=wrapCol))
        lines <- apply(base::matrix(c(HMS,msg), ncol=2),
                       1,function(x) { (sprintf("%s | %s", x[1],x[2])) })
        lines <- c(lines, sprintf("  Total: %s", tidyTime(tot)))
        paste(c(head, lines, ""), collapse="\n")
    },

    colNameToFunc = function( help=FALSE ) {
        "Internal utility, generates list-of-lists that maps color names to crayon functions"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                     names(.refClassDef@contains)) )
        use <- if(is.null(EvLogObj)) { .self } else { EvLogObj }
        if (!is.null(use$colMap) && !is.null(use$colMap$FG)) return (use$colMap)
        ## Was difficult to juggle referencing colors by function name
        ## when you can't be sure the user has installed
        ## crayon. Instead, make a named lookup of crayon functions,
        ## which will then be used by $colorize() to get() the correct
        ## function, provided it exists().
        myNames <- c("black", "red", "green", "yellow", "blue",
                     "magenta", "cyan", "white", "silver",
                     "gray", "purple", "lightblue")
        fgNames <- myNames
        names(fgNames) <- myNames
        ## I have included some aliases, remap to the R/ANSI names:
        fgNames[ "gray" ]      <- "silver"
        fgNames[ "purple" ]    <- "magenta"
        fgNames[ "lightblue" ] <- "cyan"
        ## Background color is the same, but with capitalized first
        ## letter and a "bg" prefix, eg "bgYellow":
        bgNames <- vapply(fgNames, function (x) {
            paste("bg", toupper(substr(x,1,1)),
                  substr(x,2,nchar(x)), sep="") }, "")
        bgNames[ "silver" ] <- "bgWhite"
        bgNames[ "gray" ]   <- "bgWhite"
        ## Initially I tried to manage this as a simple string->string
        ## lookup, and then converted each string to a function
        ## on-the-fly with get(). However, get() will need package
        ## crayon to be in the search space, and when EventLogger is
        ## used via inheritance that seems to not be the case. So I am
        ## going to instead use this as a string->function lookup, and
        ## evaluate each function here.
        cm <- list(FG=fgNames, BG=bgNames)
        cf <- list()
        for (typ in names(cm)) {
            cf[[typ]] <- list()
            for (nm in names(cm[[typ]])) {
                cc <- cm[[typ]][ nm ]
                ## eval() in R: https://stackoverflow.com/a/1743796
                cf[[typ]][[tolower(nm)]] <-
                    eval(parse(text=paste("crayon::", cc, sep="")))
            }
        }
        use$colMap <- cf
    },

    fieldDescriptions = function(help=FALSE) {
        "A static list of brief descriptions for each field in this object"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                              names(.refClassDef@contains)) )
        list("useCol"   = "Flag indicating if messages should be colorized",
             "log"      = "data.table storing event messages and times",
             "vb"       = "Flag indicating if messages should be shown",
             "colMap"   = "List mapping color names to crayon methods",
             "EvLogObj" = "Optional, another EventLogger object (for cross-object sharing)",
             "varName"  = "Extracted variable name associated with this object")

    },

    annotateFields = function(help=FALSE) {
        "Update object fields to include attributes with brief descriptions"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                              names(.refClassDef@contains)) )
        ## If the descriptions are not defined then nothing to be done:
        if (!is.function(.self[["fieldDescriptions"]])) return(FALSE)
        myClassName <- class(.self)
        fields <- fieldDescriptions()
        hfmt <- paste0(" help('%s', '", myClassName,
                       "') # More information on field ")
        for (fld in names(fields)) {
            if (is.null(.self[[fld]])) next # Can't attribute NULL
            ## The [[ accessor seems to work for fields?
            attr(.self[[fld]], "Description") <- fields[[fld]]
            attr(.self[[fld]], "Help") <- sprintf(hfmt,fld)
        }
        TRUE
    },

    help = function (color=NULL, help=FALSE) {
        "Display high-level help about all object methods"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                              names(.refClassDef@contains)) )
        sections <- list(
            "Logging" = c("showLog", "message", "dateMessage", "actionMessage",
            "debugMessage", "err", "logText"),
            "Formatting" = c("colorize", "useColor", "colorMap", "tidyTime",
                             "verbose"),
            "Utility Methods" = c("colNameToFunc", "fieldDescriptions")
            )
        .showHelp(sections, 'myEvtLogger', color=color)
    },

    .showHelp = function (sections=list(), genericName='myObject', color=NULL,
    help=FALSE) {
        "Construct text for the help() method"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self),
                                              names(.refClassDef@contains)) )
        if (is.null(color)) color <- useColor() # Use EventLogger setting
        doCol   <- if (color) { .self$colorize } else { function(x, ...) x }
        objName <- .self$.selfVarName(genericName)
        whtName <- doCol(objName, "white")
        comCol  <- "yellow"

        ## Figure out all available methods:
        myClassName <- class(.self)
        myClass     <- methods::getRefClass(myClassName)
        allMeth     <- myClass$methods()
        ## Subtract out the generic ones:
        allMeth <- setdiff(allMeth, c("callSuper", "copy", "export", "field", "getClass", "getRefClass", "import", "initFields", ".objectPackage", ".objectParent", "show", "trace", "untrace", "usingMethods"))
        ## Subtract out some superclasses
        allMeth <- allMeth[ !grepl('#', allMeth) ]

        ## Junk drawer section for methods not defined in `sections`:
        sections[["Other Methods"]] <- setdiff(allMeth, unname(unlist(sections)))

        ## Basic header:
        txt <- sprintf("
%s
?%s %s
%s                 %s
%s$showLog()       %s
str(%s, max.lev=3) %s
",
doCol(sprintf("###
### %s Help - call the below commands for more details
###", myClassName),"magenta"),
myClassName,
doCol("# Built-in documentation on the class", comCol),
whtName, doCol("# Summary report of the object", comCol),
whtName, doCol("# Show logged events with timing", comCol),
whtName, doCol("# Inspect the object structure", comCol))

        noHelp <- c()
        ## Add snippets for each method, broken down by section
        for (sec in names(sections)) {
            if (sec == "SKIP") next
            txt <- c(txt, doCol(paste("\n############\n###", sec, "\n"),comCol))
            meths <- sections[[ sec ]]
            for (meth in meths) {
                ## Going to see if we can extract the ROxygen
                ## description string from the method
                code <- utils::capture.output(myClass$methods(meth))
                ## Should not happen, but be safe:
                if (is.null(code)) next
                isHelped <- FALSE
                com <- NA
                for (line in code) {
                    ## See if 'help=FALSE' is set - indicates I have
                    ## tied it into my internalized help framework:
                    if (grepl('help\\s*=\\s*FALSE', line))  isHelped <- TRUE
                    cm <- CatMisc::parenRegExp('^\\s+"(.+)"\\s*$', line)
                    if (!is.na(cm[1])) {
                        ## Found a single line quoted string, presume
                        ## it is description and stop scanning code:
                        com <- cm[1]
                        break
                    }
                }
                if (!isHelped) {
                    noHelp <- c(noHelp, meth)
                    next
                }
                txt <- c(txt, sprintf("%s$%s( help=TRUE )", whtName, meth))
                if (!is.na(com)) txt <-
                     c(txt, doCol(paste("\n    #", com), comCol))
                txt <- c(txt, "\n")
            }
        }
        noHelp <- setdiff(noHelp, 'help')
        if (length(noHelp) > 0) txt <- c(txt,
              doCol("\n### Methods lacking help\n", comCol),
              sprintf("# %s\n", strwrap(paste(noHelp, collapse=' '))))

        txt <- c(txt, doCol("\n### Object fields\n", comCol))
        ## Some fields are expected to be small structures, don't need str():
        simpleFields <- unlist(strsplit("file fromRDS filterLog setFilters lvlVal colDefs", " ")) # Just a lazy way to manage the list
        strFmt <- "str(%s$%s) %s\n"
        simFmt <- "%s$%s %s\n"
        if (is.function(.self[["fieldDescriptions"]])) {
            fields <- fieldDescriptions( )
            for (field in names(fields)) {
                fmt <- ifelse(is.element(field, simpleFields), simFmt, strFmt)
                txt <- c(txt, sprintf(fmt, whtName, field,
                                      doCol(paste("#",fields[[field]]), comCol)))
            }
            annotateFields()
        } else {
            txt <- c(txt, doCol("# No field descriptions provided", "yellow"))
        }
        base::message(paste(txt, collapse='', sep=''))
        invisible(NULL)
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

