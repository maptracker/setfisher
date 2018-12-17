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
#' @seealso \link{message}, \link{showLog}
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
#' }
#'
#' This method is used to both display a message to the terminal and
#' to record it in the \link{log}
#'
#' Several wrapper functions exist with pre-configured display options:
#'
#' \itemize{
#'   \item \link{dateMessage} Will display a datestamp
#'   \item \link{actionMessage} Prefix with '[+]', red color
#'   \item \link{debugMessage} Prefix with '[DEBUG]', white FG, blue BG
#'   \item \link{err} Prefix with 'ERROR', red FG, yellow BG
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
#' \link{dateMessage}   (just sets datestamp=TRUE)
#' \link{actionMessage} (just sets prefix='[+]' and color='red')
#' \link{debugMessage}  (sets prefix='[DEBUG]', color="white", bgcolor="blue")
#' \link{err}           (sets prefix='[ERROR]', color="red", bgcolor="yellow")
#'
#' @examples
#'
#' el <- EventLogger()
#' el$message("A generic message")
#' 
#' # Show the log, nicely formatted:
#' el
NULL

#' EventLogger Date Message
#'
#' EventLogger object method to present and record a message with a date
#'
#' @name dateMessage
#' @method dateMessage EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$dateMessage( help=TRUE )
#' 
#' myObject$dateMessage(msg="No message provided!", ...)
#'
#' }
#'
#' This is just a wrapper for \link{message} where \code{datestamp=TRUE}
#'
#' @param msg The text to display and show
#' @param \dots Passed to \link{message}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The \link{log} table, invisibly
#'
#' @seealso \link{log}, \link{message}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$dateMessage("Something auspicious just happened")
#' 
#' # Show the log, nicely formatted:
#' el
NULL

#' EventLogger Action Message
#'
#' EventLogger object method to present and record an emphatic message
#'
#' @name actionMessage
#' @method actionMessage EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$actionMessage( help=TRUE )
#' 
#' myObject$actionMessage(msg="No message provided!", prefix='[+]', color="red")
#'
#' }
#'
#' This is just a wrapper for \link{message} where \code{prefix='[+]'}
#' and \code{color="red"}
#'
#' @param msg The text to display and show
#' @param prefix Default '[+]', the text to show before the message
#' @param color Default 'red', the color of the displayed text
#' @param \dots Passed to \link{message}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The \link{log} table, invisibly
#'
#' @seealso \link{log}, \link{message}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$actionMessage("Outer airlock door not responding")
#' 
#' # Show the log, nicely formatted:
#' el
NULL

#' EventLogger Debug Message
#'
#' EventLogger object method to present and record a message for debugging
#'
#' @name debugMessage
#' @method debugMessage EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$debugMessage( help=TRUE )
#' 
#' myObject$debugMessage(msg="No message provided!", prefix='[DEBUG]',
#'                       color="white", bgcolor="blue")
#'
#' }
#'
#' This is just a wrapper for \link{message} where
#' \code{prefix='[DEBUG]'}, \code{color="white"} and
#' \code{bgcolor="blue"}
#'
#' @param msg The text to display and show
#' @param prefix Default '[DEBUG]', the text to show before the message
#' @param color Default 'white', the color of the displayed text
#' @param bgcolor Default 'blue', the background color of the
#' displayed text
#' @param \dots Passed to \link{message}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The \link{log} table, invisibly
#'
#' @seealso \link{log}, \link{message}
#'
#' @examples
#'
#' el <- EventLogger()
#' for (cv in 1:9) {
#'    el$debugMessage(c("Chevron", cv, "locked ..."))
#' }
#' 
#' # Show the log, nicely formatted:
#' el
NULL

#' EventLogger Error Message
#'
#' EventLogger object method to present and record an error message
#'
#' @name err
#' @method err EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$err( help=TRUE )
#' 
#' myObject$err(msg="No message provided!", prefix='[ERROR]',
#'              color="red", bgcolor="yellow",...)
#'
#' }
#'
#' This is just a wrapper for \link{message} where
#' \code{prefix='[ERROR]'}, \code{color="red"} and
#' \code{bgcolor="yellow"}
#'
#' @param msg The text to display and show
#' @param prefix Default '[ERROR]', the text to show before the message
#' @param color Default 'red', the color of the displayed text
#' @param bgcolor Default 'yellow', the background color of the
#' @param \dots Passed to \link{message}
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The \link{log} table, invisibly
#'
#' @seealso \link{log}, \link{message}
#'
#' @examples
#'
#' el <- EventLogger()
#' el$err("Weight is not defined")
#' # Show the log, nicely formatted:
#' el
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

#' Field Descriptions
#'
#' A list of brief descriptions for each object field
#'
#' @name fieldDescriptions
#' @method fieldDescriptions EventLogger
#'
#' @details
#'
#' \preformatted{
#' ## Method Usage:
#' myObject$fieldDescriptions( help=TRUE )
#' 
#' myObject$fieldDescriptions( )
#' }
#'
#' This method returns a simple list of descriptive text for each
#' object field. It is designed to help the user understand the role
#' of each field.
#'
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A list of character strings
#'
#' @seealso \link{help}
#'
#' @examples
#'
#' myEL <- EventLogger()
#' myEL$fieldDescriptions()
#' 
NULL

