#' SetFisher Logger
#'
#' Utility class for messaging and recording of events
#'
#' @field log The data.table holding log messages
#' @field useColor Logical flag to indicate if color should be used in messaging
#' @field verbose Logical flag indicating if verbose messaging should be active
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom methods setRefClass new
#' @import crayon
#' 
#' @export SetFisherLogger
#' @exportClass SetFisherLogger
#'

SetFisherLogger <-
    setRefClass("SetFisherLogger",
                fields = list(
                    useColor = "logical",
                    log      = "data.table",
                    verbose  = "logical"
                    ),
                )

setFisherColorMap <- NULL

SetFisherLogger$methods(
    
    initialize = function(...,
        useColor = TRUE,
        verbose  = TRUE,
        log      = data.table( Date = Sys.time(),
            Message = "Log initialized", key = "Date") ) {
        callSuper(..., useColor = useColor, verbose = verbose, log = log)
       .setSetFisherColorMap( useColor )
    },

    colorize = function(msg = "", color = NULL, bgcolor = NULL) {
        ## Since setFisherColorMap will not be set up when useColor is false,
        ## this should return uncolorized text when appropriate
        if (!is.character(msg)) msg <- as.character(msg)
        if (is.def(color)   && !is.function(color))
            color <- setFisherColorMap$FG[[ tolower(color) ]]
        if (is.def(color) && exists(color)) msg <- get(color)(msg)
        if (is.def(bgcolor) && !is.function(bgcolor))
            bgcolor <- setFisherColorMap$BG[[ tolower(bgcolor) ]]
        if (is.def(bgcolor)) msg <- get(bgcolor)(msg)
        msg
        ## Getting a variable by variable string name: get()
        ## https://stackoverflow.com/a/3971855
    },
    
    message = function(msg = "No message provided!", prefix = NULL,
        color = NULL, bgcolor = NULL, datestamp = FALSE, fatal = FALSE,
        collapse = " ") {
        ## Optional prefix
        if (!is.null(prefix)) prefix <- c(prefix, " ")
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

    dateMessage = function ( msg = "No message provided!", ... ) {
        message( ... , msg = msg, datestamp = TRUE)
    },

    actionMessage = function (msg = "No message provided!!", prefix = '[+]',
        color = "red", ...) {
        message(msg = msg, prefix = prefix, color = color, ...)
    },

    debugMessage = function (msg = "No message provided!!", prefix = '[DEBUG]',
        color = "white", bgcolor = "blue", ...) {
        message(msg = msg, prefix = prefix, color = color, ...)
    },

    err = function (msg = "No message provided!!", prefix = '[ERROR]', ...) {
        message(msg = msg, prefix = prefix, collapse = "\n",
                color = "red", bgcolor = "yellow", ...)
    },

    tidyTime = function (x = NULL, pad = 0) {
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
    }
)

.setSetFisherColorMap <- function( useColor ) {
    if (!useColor) {
        setFisherColorMap <<- list()
    } else if (is.null(setFisherColorMap) ||
               is.null(setFisherColorMap$FG)) {
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
            setFisherColorMap <<- list(FG = fgNames, BG = bgNames)
        } else{
            ## crayon is not available
            message(c("Could not activate logging for messages:\n",
                      "  crayon is not installed. To install:\n",
                      "  install.packages('crayon')"))
            setFisherColorMap <<- list( FG = character(), BG = character() )
        }
    }
}

