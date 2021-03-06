% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethodsDocumentation.R
\name{message}
\alias{message}
\title{EventLogger Message}
\arguments{
\item{msg}{The text to display and show}

\item{prefix}{Default NULL, Optional text to display in front of
message. Will not be recorded in the log.}

\item{color}{Default NULL, foreground (text) color of message, not
logged.}

\item{bgcolor}{Default NUll, background color of message, not
logged.}

\item{datestamp}{Default FALSE; If TRUE, then a datestamp will be
displayed as well. Datestamps are always recorded in the
\link{log}, regardless of this value.}

\item{fatal}{Default FALSE; If TRUE, then stop() execution as well.}

\item{collapse}{Default '', text to use when collapsing msg vector.}

\item{help}{Default FALSE. If TRUE, show this help and perform no
other actions.}
}
\value{
The \link{log} table, invisibly
}
\description{
EventLogger object method to present and record a message
}
\details{
\preformatted{
## Method Usage:
myObject$message( help=TRUE )

myObject$message(msg="No message provided!", prefix=NULL,
                 color=NULL, bgcolor=NULL, datestamp=FALSE,
                 fatal=FALSE, collapse=" ")

}

This method is used to both display a message to the terminal and
to record it in the \link{log}

Several wrapper functions exist with pre-configured display options:

\itemize{
  \item \link{dateMessage} Will display a datestamp
  \item \link{actionMessage} Prefix with '[+]', red color
  \item \link{debugMessage} Prefix with '[DEBUG]', white FG, blue BG
  \item \link{err} Prefix with 'ERROR', red FG, yellow BG
}

Bear in mind that for all functions the appearance is merely
cosmetic - only the contents of \code{msg} will go into the log.
}
\examples{

el <- EventLogger()
el$message("A generic message")

# Show the log, nicely formatted:
el
}
\seealso{
\link{log}
\link{dateMessage}   (just sets datestamp=TRUE)
\link{actionMessage} (just sets prefix='[+]' and color='red')
\link{debugMessage}  (sets prefix='[DEBUG]', color="white", bgcolor="blue")
\link{err}           (sets prefix='[ERROR]', color="red", bgcolor="yellow")
}
