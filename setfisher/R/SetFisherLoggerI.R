SetFisherLoggerI <-
    setRefClass("SetFisherLoggerI")

## Wrapper calls to child objects. Seems like there must be a more
## elegant way to do this...
SetFisherLoggerI$methods(
    message       = function (...) { log$message(...)       },
    actionMessage = function (...) { log$actionMessage(...) },
    dateMessage   = function (...) { log$dateMessage(...)   },
    debugMessage  = function (...) { log$debugMessage(...)  },
    err           = function (...) { log$err(...)           },
    colorize      = function (...) { log$colorize(...)      },
    tidyTime      = function (...) { log$tidyTime(...)      }
)

