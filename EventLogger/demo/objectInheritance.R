## Example of how EventLogger can be inheritted by your own
## ReferenceClass objects, and how to share the log table between
## multiple different objects

### CAVEAT EMPTOR - I am not convinced that I am combining
### ReferenceClasses with Inheritance properly. They do *not* mix well
### together. In particular, the inheriting object does *NOT* import
### packages imported by the inherited class.

## Create a simple RefClass object that inherits (contains) EventLogger:
foo <- setRefClass("foo",
  fields   = list( x = 'numeric' ),
  contains = c("EventLogger"))

## Set the foo-Class methods:
foo$methods(
    initialize = function( x=1, ... ) {
        callSuper(x=x, ...)
    },
    set_x = function( val ) {
        x <<- val
        actionMessage(c("Set x:", val))
    },
    del_x = function() {
        x <<- as.numeric(NA)
        ## Note that this is NOT base::message(), but the RC object version:
        message("Cleared x", color='magenta', prefix='[-]')
    })

## Create a new method and manipulate it
z <- foo()
z$set_x(10)
z$set_x(3.14)
z$del_x()

## Create a second object that will share the same log as z:
q <- foo(log=z)
q$actionMessage("Hello from Q")
q$set_x(1234)
q$message(c("z =", z$x, "and q =", q$x))

## Show the log, using z
z$log

## q has a $log object, but it will alert you that it is not where
## events are being written:
q$log

## Instead, use the $EvLogObj field (which in this case is z):
q$EvLogObj$log
