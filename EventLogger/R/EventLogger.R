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
#' @field vb Logical flag indicating if verbose messaging should be
#'     active
#' @field EvLogObj An optional external EventLogger object. Used to
#' centralize log information across multiple inheriting objects in
#' one place.
#'
#' @importFrom methods new setRefClass
#' @importFrom CatMisc is.def methodHelp
#' @importFrom data.table data.table rbindlist
#' @importClassesFrom CatMisc RefClassHelper
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
                log      = "data.table",
                vb       = "logical",
                EvLogObj = "ANY"),
                contains = c("RefClassHelper") 
                )


EventLogger$methods(
    
    initialize = function(verbose=TRUE, log=NULL, ...) {
        "Create a new object using EventLogger()"
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
        ## Manage verbose flag through methods
        .self$verbose(verbose)
        callSuper(...)
    },

    help = function (color=NULL, help=FALSE) {
        "Display high-level help about all object methods"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        sections <- list(
            "Logging / Messaging" = c("message", "dateMessage",
            "actionMessage", "debugMessage", "err",
            "verbose", "logText", "showLog"),
            "Utility Methods" = c("tidyTime")
            )
        showHelp(sections, 'RefClassHelper', color=color)
    },

    fieldDescriptions = function(help=FALSE) {
        "A static list of brief descriptions for each field in this object"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        list(
            "log"      = "data.table containing actual logged data",
            "vb"       = "Verbosity flag, set with $verbose()", 
            "EvLogObj" = "Optional 'external' EventLogger object, for log centralization")
    },

    show = function (...) {
        "A wrapper for showLog"
        showLog( ... )
    },

    message = function(msg="No message provided!", prefix=NULL,
                       color=NULL, bgcolor=NULL, datestamp=FALSE,
                       fatal=FALSE, collapse=" ", help=FALSE) {
        
        "Display a formatted message, and store it in the log table"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
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
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        message(msg=msg, datestamp=TRUE, ...)
    },

    actionMessage = function (msg="No message provided!", prefix='[+]',
        color="red", help=FALSE, ...) {
        "Calls message with a '[+]' prefix and red coloring"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        message(msg=msg, prefix=prefix, color=color, ...)
    },

    debugMessage = function (msg="No message provided!", prefix='[DEBUG]',
        color="white", bgcolor="blue", help=FALSE, ...) {
        "Calls message with a '[DEBUG]' prefix and white/blue coloring"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        message(msg=msg, prefix=prefix, color=color, bgcolor=bgcolor, ...)
    },

    err = function (msg="No message provided!", prefix='[ERROR]',
    color="red", bgcolor="yellow", help=FALSE, ...) {
        "Calls message with an '[ERROR]' prefix and red/yellow coloring"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        message(msg=msg, prefix=prefix, collapse="\n",
                color=color, bgcolor=bgcolor, ...)
    },

    verbose = function(newval=NULL, help=FALSE) {
        "Get or set the flag determining if messages are displayed"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
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
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
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
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
        cat( .self$logText(...) )
    },

    logText = function (width=0.7 * getOption("width"),
                        relative=TRUE, pad=11, n=0, help=FALSE) {
        "Formats the log as a human readable two-column table"
        if (help) return( CatMisc::methodHelp(match.call(), class(.self) ) )
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
    }

)

