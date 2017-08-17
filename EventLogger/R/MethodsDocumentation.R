## ROxygen provides only rudimentary support for documenting Reference
## Class object methods (as of 2017). It will recognized an unassigned
## character string at the beginning of the method function as a
## description, and will parse the parameters. This class is complex
## enough that I wanted more formalized documentation. The blocks
## below define help topics for both RefClass fields and methods:

## * The method name is being put in @name and/or @aliases. This will
##   likely be interpreted by the help system as a simple function in the
##   namespace, but I don't think that will cause problems. It allows
##   \link{}s to be built between topics.
## * I am also specifying the methods under the @method section. This
##   will likely be perceived as an S3/S4 function in the class. I don't
##   _think_ that's an issue, but it might be.
## * It's effectively impossible to use the @usage section, since R_CMD
##   check becomes very unhappy with attempts to formalize usage in
##   'true' object fasion. Instead, usage is manually smuggled into the
##   @details section in a \preformatted{} block.

#' EventLogger Use Color Flag
#'
#' Internal EventLogger field indicating if output should be colorized
#'
#' @name useCol
#'
#' @details
#'
#' A logical flag, if TRUE it indicates that messaeges should be
#' colorized with the \link[crayon]{crayon} package.
#'
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' ## Instead, use the \link{useColor} method to check and alter the value
#' }
#'
#' @return A single logical value
#'
#' @seealso \link{useColor}, \link{message}
#'
#' @examples
#' 
#' el <- EventLogger( )
#' el$err("Something is wrong!")
#' el$useColor(FALSE)
#' el$err("It' still wrong!")
NULL

#' EventLogger Log
#'
#' Internal EventLogger field holding the data.table of events
#'
#' @name log
#'
#' @details
#'
#' A \link{data.table} (data.frame-compliant) table with two columns:
#' \code{$Date} (which is automatically populated at time of event)
#' and \code{$Message}
#'
#' This field is intended to be accessed, but alterations (additions)
#' should be performed using the \link{message} functions
#' 
#' @return A data.table
#'
#' @seealso \link{message}
#'
#' @examples
#' 
#' el <- EventLogger( )
#' el$message("An entry")
#' el$message("A second entry")
#' el$log
NULL

#' EventLogger Verbosity Flag
#'
#' Internal EventLogger field holding the data.table of events
#'
#' @name vb
#'
#' @details
#'
#' A logical flag setting if EventLogger is verbose or not.
#' 
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' ## Instead, use the \link{verbose} method to check and alter the value
#' }
#' 
#' @return A logical value
#'
#' @seealso \link{message}
#'
#' @examples
#' 
#' el <- EventLogger( )
#' el$actionMessage("HELLO")
#' el$verbose( FALSE )
#' el$message("please no shouting")
#' el
NULL

#' EventLogger Color Map
#'
#' Internal EventLogger field to manage colorizing parameters
#'
#' @name colMap
#'
#' @details
#'
#' This is a structured list holding color names and
#' \link[crayon]{crayon} functions. It helps manage selection of
#' relevant colorizing functions from a set of prameters.
#' 
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' ## ALTERING IT MAY RESULT IN CODE INSTABILITY
#' }
#' 
#' @return A list of lists
#'
#' @seealso \link{colorMap}, \link{colNameToFunc}, \link{colorize}
NULL

#' EventLogger Object
#'
#' Internal EventLogger field pointing to another object that holds a
#' shared \link{log}
#'
#' @name EvLogObj
#'
#' @details
#'
#' This package was designed to be inherited by other ReferenceClass
#' objects, and was also designed such that many objects could share a
#' common log. An object that inherits (contains) EventLogger can
#' choose to reference a different object with this field. If so,
#' operations that would normally alter the EventLogger fields will
#' instead act on the referenced object.
#'
#' The other object is defined by the \code{log} parameter on creation.
#' 
#' \preformatted{
#' ## NORMALLY YOU WILL NOT WANT TO ACCESS THIS FIELD DIRECTLY
#' ## ALTERING IT MAY RESULT IN CODE INSTABILITY
#' }
#' 
#' @return An empty field if not set, or an EventLogger-compliant object
NULL


#' EventLogger Message
#'
#' EventLogger object method to present and record a message
#'
#' @name message
#' @aliases dateMessage actionMessage debugMessage err
#' @method message EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$message( help=TRUE )
#' 
#' myObject$message(msg="No message provided!", prefix=NULL,
#'                  color=NULL, bgcolor=NULL, datestamp=FALSE,
#'                  fatal=FALSE, collapse=" ")
#'
#' myObject$dateMessage( msg )
#' myObject$actionMessage( msg )
#' myObject$debugMessage(msg)
#' myObject$err(msg)
#' }
#'
#' This method is used to both display a message to the terminal and
#' to record it in the \link{log}
#'
#' Several wrapper functions exist with pre-configured display options:
#'
#' \itemize{
#'   \item dateMessage Will display a datestamp
#'   \item actionMessage Prefix with '[+]', red color
#'   \item debugMessage Prefix with '[DEBUG]', white FG, blue BG
#'   \item err Prefix with 'ERROR', red FG, yellow BG
#' }
#'
#' Bear in mind that for all functions the appearance is merely
#' cosmetic - only the contents of \code{msg} will go into the log.
#'
#' @param msg The text to display and show
#' @param prefix Default NULL, Optional text to display in front of
#'     message. Will not be recorded in the log.
#' @param color Default NULL, foreground (text) color of message, not
#'     logged.
#' @param bgcolor Default NUll, background color of message, not
#'     logged.
#' @param datestamp Default FALSE; If TRUE, then a datestamp will be
#'     displayed as well. Datestamps are always recorded in the
#'     \link{log}, regardless of this value.
#' @param fatal Default FALSE; If TRUE, then stop() execution as well.
#' @param collapse Default '', text to use when collapsing msg vector.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The \link{log} table, invisibly
#'
#' @seealso \link{log}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$message("A generic message")
#' 
#' # Show the log, nicely formatted:
#' el
NULL

#' Map Color
#'
#' EventLogger object method to turn a string into a crayon color function
#'
#' @name colorMap
#' @method colorMap EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$colorMap( help=TRUE )
#' 
#' myObject$colorMap( color, bg=FALSE )
#' }
#'
#' Takes a color name as input, returns a \link[crayon]{crayon}
#' function. This is primarily a utility function.
#' 
#' @param color Required, a string describing a color, eg
#'     "magenta". Can also be a function reference
#' @param bg Default FALSE. If TRUE, will pick the background color
#'     corresponding to the name.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A function from the \link[crayon]{crayon}
#'
#' @seealso \link{colorize}, \link[crayon]{crayon},
#'     \link{colNameToFunc}
#'
#' @examples
#'
#' el <- EventLogger()
#' f <- el$colorMap("red")
#' base::message(f("Quick Brown Fox"))
#' f <- el$colorMap("cyan", TRUE)
#' base::message(f("Lazy Dog"))
NULL

#' Colorize
#'
#' EventLogger object method to colorize a string
#'
#' @name colorize
#' @method colorize EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$colorize( help=TRUE )
#' 
#' myObject$colorize( msg="", color=NULL, bgcolor=NULL )
#' }
#'
#' Takes a character vector and applies foregrand and/or background
#' color to it.
#' 
#' @param msg Default "". A character vector of text to colorize
#' @param color Default NULL. Optional text (eg "yellow")
#'     corresponding to the foreground color of the text
#' @param bgcolor Default NULL. Optional text (eg "silver")
#'     corresponding to the background color of the text
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character vector with ANSI color codes injected by
#'     \link[crayon]{crayon}
#'
#' @seealso \link[crayon]{crayon}, \link{colNameToFunc},
#'     \url{https://en.wikipedia.org/wiki/ANSI_escape_code#Colors}
#'
#' @examples
#'
#' el <- EventLogger()
#' x <- el$colorize(c("This", "That"), "green")
#' x
#' message(paste(x, collapse=" ... and ... "))
NULL

#' Use Color
#'
#' EventLogger object method to get/set colorization flag
#'
#' @name useColor
#' @method useColor EventLogger
#'
#' @details
#'
#' \preformatted{#' ## Method Usage:
#' myObject$useColor( help=TRUE )
#' 
#' myObject$useColor( newval=NULL )
#' }
#'
#' Gets or sets the flag determining if messages should be colorized
#' 
#' @param newval Default NULL. If provided and can be made logical,
#'     will set the flag. Inability to cast as logical will emit a
#'     non-fatal error.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A single logical value, invisibly
#'
#' @seealso \link{message}, \link{colorize}, \link{useCol}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$actionMessage("This is exciting!")
#' el$useColor(FALSE)
#' el$actionMessage("Still exciting, but it does not show as much")
NULL

#' EventLogger Verbosity
#'
#' EventLogger object method to get/set verbosity flag
#'
#' @name verbose
#' @method verbose EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$verbose( help=TRUE )
#' 
#' myObject$verbose( newval=NULL )
#' }
#'
#' Gets or sets the flag determining if messages should be displayed
#' to the terminal, or just logged.
#' 
#' @param newval Default NULL. If provided and can be made logical,
#'     will set the flag. Inability to cast as logical will emit a
#'     non-fatal error.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A single logical value, invisibly
#'
#' @seealso \link{message}, \link{colorize}, \link{vb}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$err("Please be aware that something has gone wrong")
#' el$verbose(FALSE)
#' el$err("Another problem! But not on your screen. You'll need to check the log")
#' el
NULL

#' Tidy Time
#'
#' EventLogger object method to pretty-format a time difference
#'
#' @name tidyTime
#' @method tidyTime EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$tidyTime( help=TRUE )
#' 
#' myObject$tidyTime( (x=NULL, pad=0 )
#' }
#'
#' Takes a character vector and applies foregrand and/or background
#' color to it. This method is used by \link{showLog} to highlight
#' shorter or longer time frames.
#' 
#' @param x Default NULL. Expects a single numeric value in seconds
#' @param pad Default 0, a minimum width that the final string should
#'     occupy
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character string with ANSI color codes injected by
#'     \link[crayon]{crayon}
#'
#' @seealso \link{showLog}, \link{logText}
#'
#' @examples
#'
#' el <- EventLogger()
#' base::message( el$tidyTime(0.00013) )
#' base::message( el$tidyTime(35232) )
NULL

#' Show Log
#'
#' EventLogger object method to pretty-print the event log
#'
#' @name showLog
#' @method showLog EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$showLog( help=TRUE )
#' 
#' myObject$showLog( ... )
#' }
#'
#' A simple wrapper for logText(), the contents of which are displayed
#' in the shell using \code{cat}. This method is auto-invoked if an
#' EventLogger object is evaluated in the shell (and there are no
#' other ReferenceClass classes contained in the object that would
#' take precedence)
#' 
#' @param ... Will be passed to \link{logText}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return NULL
#'
#' @seealso \link{logText}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$message("A sample message")
#' el$showLog()
NULL

#' Log Text
#'
#' EventLogger object method to generate pretty-formatted text of the log
#'
#' @name logText
#' @method logText EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$logText( help=TRUE )
#' 
#' myObject$logText( width=0.7 * getOption("width"),
#'                   relative=TRUE, pad=11, n=0 )
#' }
#'
#' This method will parse the \link{log} data.table and generate a
#' human-friendly, colorized table of events. The left column will
#' report the time difference between log events, intending to help
#' gauge how long individual events are taking. The right column will
#' be event message text, strwrap()'ed to the user's \code{width} option.
#'
#' This method will return the 'raw' character strings of the text. In
#' general you will likely wish to instead use \link{showLog} to have
#' the message shown in a way such that the ANSI color codes are
#' properly evaluated for display.
#' 
#' @param widthDefault 70% of the \code{option('width')}. This is the
#'     length to be used when strwrap()ing the event text
#' @param relative - Default TRUE, will show the time elapsed between
#'     events. If FALSE will show absolute time stamps
#' @param pad - Default 11, character padding used to align the time
#'     column
#' @param n Default 0, if greater will only show that number of most
#'     recent events
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return NULL
#'
#' @seealso \link{showLog}, \link{message}, \link{log}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$message("A sample message")
#' el$logText()
NULL

#' Color Name to Crayon Function
#'
#' EventLogger object method to generate a map of color names to
#' crayon functions
#'
#' @name colNameToFunc
#' @method colNameToFunc EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$colNameToFunc( help=TRUE )
#' 
#' myObject$colNameToFunc( )
#' }
#'
#' This list generated by this method is simply a lookup list-of-lists
#' used to turn names ("cyan", "black") into \link[crayon]{crayon}
#' colorizing functions. It is utilized by \link{colorize} to convert
#' the user's color parameters into the appropriate colorizing
#' functions.
#'
#' The function is run only when \link{colMap} has not already been
#' set up. Afterwards the resulting list is stored in \link{colMap}
#' for cached retrieval.
#' 
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A list of lists
#'
#' @seealso \link{colMap}, \link{colorize}
NULL

