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
#' @field useColor Logical flag to indicate if color should be used in messaging
#' @field verbose Logical flag indicating if verbose messaging should be active
#' @field colMap A list generated to map color names to crayon functions
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom methods new setRefClass
#' @importClassesFrom data.table data.table
#' @import crayon
#'
#' @examples
#'
#' log <- EventLogger()
#' x <- "wobble"
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
#' log$log
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
#' z$log
#' 
#' @export EventLogger
#' @exportClass EventLogger
#'

EventLogger <-
    setRefClass("EventLogger",
                fields = list(
                    useColor = "logical",
                    log      = "data.table",
                    verbose  = "logical",
                    colMap   = "list"
                    ),
                )

EventLogger$methods(
    
    initialize = function(useColor=TRUE, verbose =TRUE, log=NULL, ...) {
        initFields(useColor=useColor, verbose=verbose, ...)
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
                    data.table(Date=Sys.time(),
                               Message="Log initialized", key="Date")
                } else { log }
        .setEventColorMap( useColor )
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
        } else if (verbose) {
            base::message(m)
        }
        ## Add entry to log table
        log <<- rbindlist(list(log, data.table(
            Date = ti, paste(msg, collapse = " "))))
        invisible(log)
    },

    colorize = function(msg = "", color = NULL, bgcolor = NULL) {
        "\\preformatted{Colorize a string with crayon. Parameters:
      msg - The string to colorize
    color - The color (eg 'red') to assign to the text
  bgcolor - The color to assign to the background
}"
        
        ## Since colMap will not be set up when useColor is false,
        ## this should return uncolorized text when appropriate
        if (!is.character(msg)) msg <- as.character(msg)
        if (.isDef(color)   && !is.function(color))
            color <- colMap$FG[[ tolower(color) ]]
        if (.isDef(color) && exists(color)) msg <- get(color)(msg)
        if (.isDef(bgcolor) && !is.function(bgcolor))
            bgcolor <- colMap$BG[[ tolower(bgcolor) ]]
        if (.isDef(bgcolor)) msg <- get(bgcolor)(msg)
        msg
        ## Getting a variable by variable string name: get()
        ## https://stackoverflow.com/a/3971855
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

    show = function ( width = 0.7 * getOption("width"),
        relative = TRUE, pad = 11, n = 0) {
        "Pretty-prints the log, including total elapsed time"
        usingMethods("tidyTime") # Needed for use in apply
        head <- colorize("Activity log:", "blue")
        ## Nicely format the log
        HMS <- log$Date
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
        msgs <- log$Message
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

    .setEventColorMap = function( useColor ) {
        if (!useColor) {
            colMap <<- list()
        } else if (is.null(colMap) || is.null(colMap$FG)) {
            ## if (require("crayon", quietly = TRUE)) {
            if (TRUE) {
                ## Was difficult to juggle referencing colors by function
                ## name when you can't be sure the user has installed
                ## crayon. Instead, make a named lookup of crayon
                ## functions, which will then be used by $colorize() to
                ## get() the correct function, provided it exists().
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
                ## letter and a "bg" prefix:
                bgNames <- vapply(fgNames, function (x) {
                    paste("bg", toupper(substr(x,1,1)),
                          substr(x,2,nchar(x)), sep="") }, "")
                colMap <<- list(FG = fgNames, BG = bgNames)
            } else{
                ## crayon is not available
                message(c("Could not activate color logging for messages:\n",
                          "  crayon is not installed. To install:\n",
                          "  install.packages('crayon')"))
                colMap <<- list( FG = character(), BG = character() )
            }
        }
        colMap
    },

    .isDef = function (x) {
        # "is defined"
        if (is.null(x)) {
            FALSE
        } else if (length(x) == 0 || all(is.na(x))) {
            FALSE
        } else {
            TRUE
        }
    }
)

