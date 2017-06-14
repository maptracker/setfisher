#' Event Logger Interface
#'
#' An interface layer designed to be inherited by a RefClass object
#' that has a 'log' field holding an EventLogger object. I am not
#' certain this is the appropriate way to expose EventLogger methods
#' to an inheriting object...
#'
#' The use case is a log table shared (hence the use of RefClass) by
#' multiple objects.
#'
#' Note that an object inheriting EventLoggerI will have internal
#' calls to \code{message()} using \code{EventLogger::message()}, not
#' \code{base::message()}.
#'
#' @seealso \link{EventLogger} for the primary object class.
#'
#' @examples
#'
#' ## Create a simple RefClass object that inherits (contains) EventLogger:
#' foo <- setRefClass("foo",
#'   fields = list( x = 'numeric', log = 'EventLogger'),
#'   contains = c("EventLoggerI"))
#' 
#' ## Set the foo-Class methods:
#' foo$methods(
#'     initialize = function(..., x=1, log=EventLogger()) {
#'         callSuper(..., x=x, log=log)
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
#' @export EventLoggerI
#' @exportClass EventLoggerI
#'

EventLoggerI <-
    setRefClass("EventLoggerI")

EventLoggerI$methods(
    message       = function (...) { log$message(...)       },
    actionMessage = function (...) { log$actionMessage(...) },
    dateMessage   = function (...) { log$dateMessage(...)   },
    debugMessage  = function (...) { log$debugMessage(...)  },
    err           = function (...) { log$err(...)           },
    colorize      = function (...) { log$colorize(...)      },
    tidyTime      = function (...) { log$tidyTime(...)      }
)

## Surely there's a better way to do this...
