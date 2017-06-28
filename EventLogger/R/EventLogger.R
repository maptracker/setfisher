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
#' @importFrom CatMisc is.def
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
                    EvLogObj = "ANY"
                    ),
                )

EventLogger$methods(
    
    initialize = function(useColor=NULL, verbose=NULL, log=NULL, ...) {
        "\\preformatted{Create a new object using EventLogger():
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
        callSuper(...)
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

    message = function(msg = "No message provided!", prefix = NULL,
                       color = NULL, bgcolor = NULL, datestamp = FALSE, fatal = FALSE,
                       collapse = " ") {
        "\\preformatted{Display an optionally colorized message, and store
it in the log table. Parameters:
      msg - The text to display and show
   prefix - Optional text to display in front of message.
            Will not be recorded in the log.
    color - Foreground (text) color of message, not logged.
  bgcolor - Background color of message, not logged.
datestamp - If TRUE, then a datestamp will be displayed as well.
            The stamp will not be recorded with the log message; The
            $log data.table includes a 'Date' column.
    fatal - If TRUE, then stop() execution as well.
 collapse - Default '', text to use when collapsing msg vector.
}"

        if (!is.null(prefix)) prefix <- c(prefix, " ") # Optional prefix
        ti <- Sys.time()
        ds <- if (datestamp) {
            ## Append date stamp to message
            ds <- colorize(as.character(ti), "yellow")
            c(" - ", ds)
        } else {
            NULL
        }
        m <- colorize(c(prefix, paste(msg, collapse = collapse)), color)
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
                 Date = ti, paste(msg, collapse = " "))))
            invisible(EvLogObj$log)
        } else {
            log <<-
                data.table::rbindlist(list(log, data.table::data.table(
                 Date = ti, paste(msg, collapse = " "))))
            invisible(log)
        }
    },

    colorMap = function(color, bg=FALSE) {
        "\\preformatted{Picks the appropriate crayon color for a color name
    color - The name (string) of the color, *or* a function reference
       bg - Default FALSE. If TRUE, then the method will return the relevant
            background color method
}"
        
        if (!useColor() || ! CatMisc::is.def(color)) return(NA)
        if (is.function(color)) return( color )
        key    <- ifelse(bg, "BG", "FG")
        colNameToFunc()[[ key ]][[ tolower(color) ]]
    },

    colorize = function(msg = "", color = NULL, bgcolor = NULL) {
        "\\preformatted{Colorize a string with crayon. Parameters:
      msg - The string to colorize
    color - The color (eg 'red') to assign to the text
  bgcolor - The color to assign to the background
}"
        
        if (!is.character(msg)) msg <- as.character(msg)
        fgFn <- colorMap(color, FALSE)
        if (is.function(fgFn)) msg <- fgFn(msg)
        bgFn <- colorMap(bgcolor, TRUE)
        if (is.function(bgFn)) msg <- bgFn(msg)
        msg
    },
    
    useColor = function(newval=NULL) {
        "\\preformatted{Get/Set flag to use color or not. Parameters:
   newval - Optional new value. Should be logical or as.logical()-able.

     Current value is returned invisibly
}"
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

    verbose = function(newval=NULL) {
        "\\preformatted{Get/Set flag for messaging to be verbose or not. Parameters:
   newval - Optional new value. Should be logical or as.logical()-able.

     Current value is returned invisibly
}"
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

    dateMessage = function ( msg = "No message provided!", ... ) {
        "Calls message() with datestamp=TRUE"
        message(msg=msg, datestamp = TRUE, ...)
    },

    actionMessage = function (msg = "No message provided!!", prefix = '[+]',
        color = "red", ...) {
        "Calls message with a '[+]' prefix and red coloring"
        message(msg=msg, prefix=prefix, color=color, ...)
    },

    debugMessage = function (msg = "No message provided!!", prefix = '[DEBUG]',
        color = "white", bgcolor = "blue", ...) {
        "Calls message with a '[DEBUG]' prefix and white/blue coloring"
        message(msg=msg, prefix=prefix, color=color, bgcolor=bgcolor, ...)
    },

    err = function (msg = "No message provided!!", prefix = '[ERROR]', ...) {
        "Calls message with an '[ERROR]' prefix and red/yellow coloring"
        message(msg=msg, prefix=prefix, collapse="\n",
                color="red", bgcolor="yellow", ...)
    },

    tidyTime = function (x = NULL, pad = 0) {
        "\\preformatted{Reports a time interval with unit management and colorization based on overall elapsed time
        x - The time unit to be tidied, in seconds
      pad - Default 0, a minimum width that the final string should occupy
}"
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

    show = function (...) {
        "\\preformatted{Pretty-prints the log, including total elapsed time
A simple wrapper for logText(), and is auto-invoked if an EventLogger
object is evaluated in the shell
}"
        cat( .self$logText(...) )
    },

    logText = function ( width = 0.7 * getOption("width"),
                        relative = TRUE, pad = 11, n = 0) {
        "\\preformatted{Format the log data to show events and elapsed times
    width - Default 70% of the 'widht' option. The length to be used when
            strwrap()ing the event text
 relative - Default TRUE, will show the time elapsed between events. If
            FALSE will show absolute time stamps
      pad - Default 11, character padding
        n - Default 0, if greater will only show that number of most
            recent events
}"
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
                ## The last time has no delta, add an empty string
                HMS <- c( unlist(lapply(HMS, tidyTime, pad = pad)),
                         strrep(" ",pad))
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
        wrapCol = paste("\n", strrep(" ", pad + 3), sep = "")
        msg <- unlist(lapply(apply(as.matrix(msgs),
                                   1, strwrap, width = width),
                             paste, collapse = wrapCol))
        lines <- apply(matrix(c(HMS,msg), ncol = 2),
                       1,function(x) { (sprintf("%s | %s", x[1],x[2])) })
        lines <- c(lines, sprintf("  Total: %s", tidyTime(tot)))
        paste(head, lines, sep = "\n")
    },

    colNameToFunc = function( ) {
        "Internal utility, generates list-of-lists that maps color names to crayon functions"
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
        bgNames[ "gray" ] <- "bgWhite"
        ## Initially I tried to manage this as a simple string->string
        ## lookup, and then converted each string to a function
        ## on-the-fly with get(). However, get() will need package
        ## crayon to be in the search space, and when EventLogger is
        ## used via inheritance that seems to not be the case. So I am
        ## going to instead use this as a string->function lookup, and
        ## evaluate each function here.
        cm <- list(FG = fgNames, BG = bgNames)
        cf <- list()
        for (typ in names(cm)) {
            cf[[typ]] <- list()
            for (nm in cm[[typ]]) {
                ## eval() in R: https://stackoverflow.com/a/1743796
                cf[[typ]][[nm]] <-
                    eval(parse(text=paste("crayon::", nm, sep="")))
            }
        }
        use$colMap <- cf
    }
)

