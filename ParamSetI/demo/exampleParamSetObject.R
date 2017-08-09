## Example of how ParamSetI can be inheritted by your own
## ReferenceClass objects

### CAVEAT EMPTOR - I am not convinced that I am combining
### ReferenceClasses with Inheritance properly. They do *not* mix well
### together. In particular, the inheriting object does *NOT* import
### packages imported by the inherited class. From my reading of the
### Reference Class documentation, this seems to be a known limitation
### of the class structure.

## Create a simple RefClass object that inherits (contains) ParamSetI:
ParamSetExample <- setRefClass("ParamSetExample",
  fields   = list( x = 'integer' ),
  contains = c("ParamSetI"))

## Set the ParamSetExample-Class methods:
ParamSetExample$methods(
    initialize = function( x=1L, params=NULL, ... ) {
        x <<- x
        ## Set some default parameters:
        setParamList(list(inc=2L,dec=3L))
        ## Then set any user values based on params
        callSuper(params=params, ...)
        defineParameters("
inc    [integer] Increment amount
dec    [integer] Decrement amount
Prod   [numeric] A mysterious value that is not initially set
")
    },
    inc = function( val ) {
        "Increment x by the 'inc' parameter."
        v <-  param("inc")
        x <<- x + v
        debugMessage("Incremented x by",v,": now =", x)
    },
    dec = function( val ) {
        "Decrement x by the 'dec' parameter."
        v <-  param("dec")
        x <<- x - v
        debugMessage("Decremented x by",v,": now =", x)
    },
    show = function () {
        "Pretty-print the object"
        cat(sprintf("ParamSetExample object:
  x = %d
      inc: %d
      dec: %d
      Available Parameters:
", x, param("inc"), param("dec") ))
        showParameters()
    }
)

## Object is defined and ready for use
