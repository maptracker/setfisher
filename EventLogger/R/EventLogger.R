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
#' @field log A list object holding the data.table which records the
#'     log messages
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom methods new setRefClass
#' @importFrom CatMisc is.def
#' @importClassesFrom data.table data.table
#' @import crayon
#'
#' @examples
#'
#' log <- EventLogger()
#' log$message("Did something important")
#' Sys.sleep(3)
#' log$actionMessage("Something emphatic has happened")
#' Sys.sleep(1)
#' log$dateMessage("Here's a date stamp")
#' log$debugMessage("Remember to comment this out in production")
#'
#' # Pretty print the log, including an elapsed time:
#' log
#' # Expose the underlying data.table:
#' log$log$log # Everyone loves it
#'
#' # Have this class inherited by another RefClass object:
#' 
#' ## Create a simple RefClass object that inherits (contains) EventLogger:
#' foo <- setRefClass("foo",
#'   fields = list( x = 'numeric' ), contains = c("EventLogger"))
#' 
#' ## Set the foo-Class methods:
#' foo$methods(
#'     initialize = function( x=1, ... ) {
#'         initFields(x=x)
#'         callSuper(...)
#'     },
#'     set_x = function( val ) {
#'         x <<- val
#'         actionMessage(c("Set x:", val))
#'     },
#'     del_x = function() {
#'         x <<- as.numeric(NA)
#'         ## Note that this is NOT base::message(), but the RC object version:
#'         message("Cleared x", color='magenta', prefix='[-]')
#'     })
#' 
#' ## Create a new method and manipulate it
#' z <- foo()
#' z$set_x(10)
#' z$set_x(3.14)
#' z$del_x()
#' 
#' ## Show the log:
#' z$log$log
#' 
#' @export EventLogger
#' @exportClass EventLogger
#'

EventLogger <-
    setRefClass("EventLogger",
                fields = list(
                    log      = "list",
                    verbose  = "logical"
                    ),
                )

EventLogger$methods(
    
    initialize = function(useColor=TRUE, verbose=TRUE, log=NULL, ...) {
        log <<- if (is.null(log)) {
                    
### This require() statement is silly when using EventLogger directly
### - the package is imported and available when creating objects with
### EventLogger(). However, when using it as an inherited class (via
### contains=) it seems to be neccesary. If I fail to include this
### here, the generation of the default log object fails with:
###    Error in callSuper(..., x = x) : could not find function "data.table"
### It seems like I am doing something wrong setting up inheritance,
### but I can't for the life of me figure out what.
                    
                    require("data.table", quietly=TRUE)
                    require("CatMisc", quietly=TRUE)
                    list(useColor=useColor,
                         verbose=verbose,
                         log=data.table(Date=Sys.time(),
                                        Message="Log initialized", key="Date"))
                } else {
                    log
                }
        .setEventColorMap( )
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
        } else if (log$verbose) {
            base::message(m)
        }
        ## Add entry to log table
        log$log <<- rbindlist(list(log$log, data.table(
            Date = ti, paste(msg, collapse = " "))))
        invisible(log$log)
    },

    useColor = function(newval=NULL) {
        "\\preformatted{Get/Set flag to use color or not. Parameters:
   newval - Optional new value. Should be logical or as.logical()-able.
}"
        if (!is.null(newval)) {
            nv <- as.logical(newval)[1]
            if (is.na(nv)) {
                err("useColor() should be provided with a boolean argument")
            } else {
                log$useColor <<- nv
            }
        }
        log$useColor
    },

    colMap = function(color, bg=FALSE) {
        "\\preformatted{Picks the appropriate crayon color for a color name
    color - The name (string) of the color, *or* a function reference
       bg - Default FALSE. If TRUE, then the method will return the relevant
            background color method
}"
        
        if (!useColor()) return(NA)
        if (is.function(color)) return( color )
        key    <- ifelse(bg, "BG", "FG")
        fnName <- log$colMap[[ key ]][ tolower(color) ]
        if (is.def(fnName) && exists(fnName)) {
            ## Getting a variable by variable string name: get()
            ## https://stackoverflow.com/a/3971855
            get(fnName)
        } else {
            NA
        }
    },

    colorize = function(msg = "", color = NULL, bgcolor = NULL) {
        "\\preformatted{Colorize a string with crayon. Parameters:
      msg - The string to colorize
    color - The color (eg 'red') to assign to the text
  bgcolor - The color to assign to the background
}"
        
        if (!is.character(msg)) msg <- as.character(msg)
        fgFn <- colMap(color, FALSE)
        if (is.function(fgFn)) msg <- fgFn(msg)
        bgFn <- colMap(bgcolor, TRUE)
        if (is.function(bgFn)) msg <- bgFn(msg)
        msg
    },
    
    dateMessage = function ( msg = "No message provided!", ... ) {
        "\\preformatted{Calls message() with datestamp=TRUE
}"
        message(msg=msg, datestamp = TRUE, ...)
    },

    actionMessage = function (msg = "No message provided!!", prefix = '[+]',
        color = "red", ...) {
        "\\preformatted{Calls message with a '[+]' prefix and red coloring
}"
        message(msg=msg, prefix=prefix, color=color, ...)
    },

    debugMessage = function (msg = "No message provided!!", prefix = '[DEBUG]',
        color = "white", bgcolor = "blue", ...) {
        "\\preformatted{Calls message with a '[DEBUG]' prefix and white/blue coloring
}"
        message(msg=msg, prefix=prefix, color=color, bgcolor=bgcolor, ...)
    },

    err = function (msg = "No message provided!!", prefix = '[ERROR]', ...) {
        "\\preformatted{Calls message with an '[ERROR]' prefix and red/yellow coloring
}"
        message(msg=msg, prefix=prefix, collapse="\n",
                color="red", bgcolor="yellow", ...)
    },

    tidyTime = function (x = NULL, pad = 0) {
        "Reports a time interval with unit management and colorization based on overall elapsed time"
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

    showLog = function ( width = 0.7 * getOption("width"),
        relative = TRUE, pad = 11, n = 0) {
        "Pretty-prints the log, including total elapsed time"
        usingMethods("tidyTime") # Needed for use in apply
        head <- colorize("Activity log:", "blue")
        ## Nicely format the log
        HMS <- log$log$Date
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
        msgs <- log$log$Message
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
        cat(head, lines, sep = "\n")
        invisible(NULL)
    },

    show = function (...) showLog(...),

    .setEventColorMap = function( useColor ) {
        if (!is.null(log$colMap) && !is.null(log$colMap$FG))
            return (log$colMap)
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
        log$colMap <<- list(FG = fgNames, BG = bgNames)
    }
)

