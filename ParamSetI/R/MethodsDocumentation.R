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


#' Paramater
#'
#' ParamSetI object method to get/set a parameter
#'
#' @name param
#' @method param ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$param( help=TRUE )
#'
#' myObject$param(key=NA, val=NA, append=FALSE, default=NA, clobber=TRUE, 
#'                check.class=NULL, is.scalar=NULL, coerce=TRUE )
#' }
#'
#' Gets or sets parameters on the object
#'
#' @param key Default NA. The name of the parameter. If NA will simply
#'     return NA. Key names are case insensitive.
#' @param val Default NA. Optional new value. If not NA, then the
#'     parameter value will be set, with behavior modified by some of
#'     the flags below
#' @param append Default FALSE, which will cause the parameter to be
#'     set to \code{val}. If TRUE, \code{val} will be appended to the
#'     end of the current vector holding the value for \code{key}.
#' @param default Default NA. Optional value to return if the value of
#'     \code{key} is 'not defined'; see \link[CatMisc]{is.def} for
#'     values that are considered 'undefined'
#' @param clobber Default TRUE, which allows an already-set value to
#'     be replaced with \code{val}. Using FALSE is primarily used for
#'     managing default settings.
#' @param check.class Default NULL, which will check the parameter
#'     definitions and use any class found there. If the value is NA
#'     or '', then there will be no class check. Otherwise, is.class()
#'     will be tested with the provided class name against the
#'     provided \code{val} to see if it matches. If not, an error will
#'     be reported and \code{key} will not be set. The value 'percent'
#'     will be interpreted as 'numeric'.
#' @param is.scalar Default NULL. If TRUE then only \code{val[1]} will
#'     be taken. If FALSE, then all elements of \code{val} to be
#'     assigned to \code{key}. If NULL, then if val is only a single
#'     element, but matches the pattern '[^]][.+]', then it will be
#'     split as a vector - see \link{selfSplittingString}
#' @param coerce Default TRUE, which will attempt to coerce \code{val}
#'     to \code{check.class} in the event that check.class is not NA
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return The value of the parameter, or NA if the key does not exist
#'     or is unset
#'
#' @seealso \link{showParameters}, \link{paramClass}, \link[CatMisc]{is.def}
#'
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' # Get a parameter
#' pse$param("inc")
#' # Set a parameter
#' pse$param("inc", 33L)
#' 
#' # Setting a parameter to an illegal class is prevented
#' pse$param("inc", "swamp monster")
#' # ... unless you take the safety off:
#' pse$param("inc", "swamp monster", check.class=FALSE)
#' # ... or specify an allowed class:
#' pse$param("inc", "swamp monster", check.class="character")
#'
#' # coercable classes are allowed:
#' pse$param("inc", 42) # Numeric, not integer
#' # ... unless you don't allow it
#' pse$param("inc", 42, coerce=FALSE)
#' 
#' # Multiple values are normally ignored
#' pse$param("inc", 4:7)
#' # ... unless you indicate otherwise:
#' pse$param("inc", 4:7, is.scalar=FALSE)
#' # ... or you indicate that you want to append new values
#' pse$param("inc", 8:10, append=TRUE) # but only the first value!
#' # ... or you indicate *both*:
#' pse$param("inc", 11:13, is.scalar=FALSE, append=TRUE)
#' # ... or pass a string recognized by selfSplittingString
#' pse$param("inc", "[..][2..22..222]")
#' 
#' # clobber can be used to prevent overwriting a previously-set value
#' pse$param("inc", 1L, clobber=FALSE)
#'
#' # You can specify a default value if none is already set
#' pse$param("color")
#' pse$param("color", default="purple")
#' 
#' 
NULL

#' All Parameters
#'
#' ParamSetI object method to list all parameter keys
#'
#' @name allParams
#' @method allParams ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$allParams( help=TRUE )
#'
#' myObject$allParams()
#' }
#'
#' Lists all "known" parameter key names. This includes parameters
#' that have been set, as well as those that are unset but have been
#' described.
#'
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A character vector of parameter key names
#'
#' @seealso \link{hasParam}, \link{showParameters}, \link{param}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' pse$allParams()
#' pse$param("puppy", "fuzzy") # Add a new parameter
#' pse$allParams()
#'
NULL

#' Has Parameter
#'
#' ParamSetI object method to test if an object has a parameter or not
#'
#' @name hasParam
#' @method hasParam ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$hasParam( help=TRUE )
#'
#' myObject$hasParam( key=NULL )
#' }
#'
#' Checks if the provided key(s) are defined in the parameter
#' set. Case/capitalization is ignored.
#'
#' @param key Default NULL. One or more parameter names
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A logical vector with TRUE indicating presence of the
#'     corresponding key in the parameter set
#'
#' @seealso \link{allParams}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' pse$showParameters()
#' pse$hasParam("inc")
#' pse$hasParam("maximum operating depth")
#'
NULL

#' Parameter Definition
#'
#' ParamSetI object method to get/set the definition (description) of a parameter
#'
#' @name paramDefinition
#' @method paramDefinition ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$paramDefinition( help=TRUE )
#'
#' myObject$paramDefinition( key, val=NULL )
#' }
#'
#' Parameters can have optional descriptions to help inform users of
#' their usage. This function returns these definitions, or allows
#' them to be set. Deffinitions will be displayed in some output
#' functions, such as \link{showParameters}.
#'
#' For setting large numbers of parameters, also consider
#' \link{defineParameters}, which can parse descriptions from a block of
#' text.
#'
#' @param key Required, one or more parameter names
#' @param val Default NULL. Optional character vector of new parameter
#'     definitions
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A named character vector of definitions, with names
#'     corresponding to keys
#'
#' @seealso \link{defineParameters}, \link{allParams}, \link{showParameters}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' # Get the definition for a parameter
#' pse$paramDefinition("inc")
#' # Set the definitions for some new parameters
#' pse$paramDefinition(c("color", "velocity"),
#'                     c("The color of the object", "The speed of the object"))
#' # Show all descriptions
#' pse$paramDefinition( pse$allParams() )
#'
NULL

#' Parameter Class
#'
#' ParamSetI object method to get/set class restrictions of a parameter
#'
#' @name paramClass
#' @method paramClass ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$paramClass( help=TRUE )
#'
#' myObject$paramClass( key, val=NULL )
#' }
#'
#' Parameters can have optional class restrictions (storage modes, eg
#' 'integer' or 'logical') assigned to them. These serve as sanity
#' checks when calling \link{param}; Unless \code{check.class} is set
#' to FALSE, \link{param} will refuse to set a mismatched parameter
#' value.
#'
#' For setting large numbers of parameters, also consider
#' \link{defineParameters}, which can parse a class restriction from a
#' block of text, or \link{setParamList}, which can set parameters
#' from a list and provide \code{check.class} as dots to \link{param}.
#'
#' @param key Required, one or more parameter names
#' @param val Default NULL. Optional character vector of new parameter
#'     class restrictions
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A named character vector of class restrictions, with names
#'     corresponding to keys
#'
#' @seealso \link{defineParameters}, \link{showParameters}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' # Get the class restriction for a parameter
#' pse$paramClass("inc")
#' 
#' # Set class restrictions for some new parameters
#' pse$paramClass(c("count", "isOnFire"),
#'                c("integer", "logical"))
#' # New class restrictions prevents illegal assignments
#' pse$param("isOnFire", "OMG YES")
#' pse$param("isOnFire", TRUE)
#'
NULL

#' Parameter Name
#'
#' ParamSetI object method to get/set the pretty-print name of a parameter
#'
#' @name paramName
#' @method paramName ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$paramName( help=TRUE )
#'
#' myObject$paramName( key, val=NULL )
#' }
#'
#' A mostly-aesthetic function to configure the case/capitalization of
#' parameter names. Internally the package ignores case assigned to
#' parameters. However, for display purposes the capitalization can be
#' specifically set. This method allows a particular capitalization to
#' be defined.
#'
#' @param key Required, one or more parameter names
#' @param val Default NULL. Optional character vector of new parameter
#'     names. \code{val} can only differ from \code{key} in the
#'     capitalization, otherwise an error will be thrown.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @return A named character vector of class restrictions, with names
#'     corresponding to keys
#'
#' @seealso \link{defineParameters}, \link{showParameters}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' # Get the pretty print name for a parameter
#' pse$paramClass("inc")
#' 
#' # That's boring. SHOUT IT
#' pse$paramClass("inc", "INC")
#' 
#' # Mostly handy for CamelCase. Note that key can be any case
#' pse$param("ADARKANDSTORMYNIGHT", "aDarkAndStormyNight")
NULL

#' Set Param List
#'
#' ParamSetI object method to set parameters in bulk using a list object
#'
#' @name setParamList
#' @method setParamList ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$setParamList( help=TRUE )
#'
#' myObject$setParamList( params=NULL, ... )
#' }
#'
#' A mechanism to bulk set parameters using a list. The list names
#' will be taken as keys, and their contents as values.
#'
#' @param params Default NULL, which will be a no-op and return
#'     NA. Otherwise, a list holding the parameters, with the names
#'     representing key names and the contents the values to assign to
#'     each key
#' @param ... Any other parameters will be passed to param() for each
#'     key/value pair
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @seealso \link{defineParameters}, \link{param},
#'     \link{showParameters}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' # What is already set:
#' pse$showParameters()
#' 
#' # Reset some existing parameters:
#' pse$setParamList( list( inc=7L, dec=2L ) )
#' # Set some new parameters:
#' pse$setParamList( list(color="green", speed="88mph") )
#' # Check current values
#' pse$showParameters()
#'
#' # Can unset clobber to set defaults without overwriting current values
#' pse$setParamList( list( inc=999L, dec=888L ), clobber=FALSE )
#' pse$showParameters()
#' 
NULL

#' Define Parameters Set Param List
#'
#' ParamSetI object method to bulk set definitions (descriptions) for parameters
#'
#' @name defineParameters
#' @method defineParameters ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$defineParameters( help=TRUE )
#'
#' myObject$defineParameters( x )
#' }
#'
#' A mechanism to bulk set parameters using a block of text. The text
#' block will be split into lines and parsed with regular expressions
#' to extract names, class restrictions and definitions. The method
#' does not define values.
#'
#' @param x Required, the block of text to parse (character)
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @seealso \link{param}, \link{setParamList}, \link{showParameters}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' # Define four new parameters, three of which have class
#' # restrictions, one does not:
#' pse$defineParameters("
#'    color [character] The color of the object, a name or hex value
#'   weight [numeric]   The mass of the object, in kilograms
#'  inStock [logical]   Flag to indicate if inventory is available
#'     misc Random information, see Jacob in logistics for more info
#' ")
#'
#' # Set a few of these
#' pse$param("weight", 74.3)
#' pse$param("instock", "yes") # Oops, violates class restriction
#' pse$param("misc", "[/][Stack on lower shelf/No hooks]")
#' 
#' # Show all the definitions
#' pse$showParameters( na.rm=FALSE )
#' 
NULL

#' Show Parameters
#'
#' ParamSetI object method to pretty-print parameters and their values
#'
#' @name showParameters
#' @method showParameters ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$showParameters( help=TRUE )
#'
#' myObject$showParameters( na.rm=TRUE )
#' }
#'
#' Show parameters in a format designed for human
#' consumption. Parameters will be shown wrapped in a \code{$param()}
#' call to illustrate how to access and change the parameter. If a
#' definition has been set, it will be shown as a comment line under
#' the parameter value.
#'
#' @param na.rm Default TRUE, which will exclude parameters that are
#'     not "defined" (see \link[CatMisc]{is.def} for values that are
#'     considered 'undefined'). Setting to FALSE shows all parameters
#'     available from \link{allParams}.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @seealso \link{param}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' pse <- ParamSetExample( params=list(inc=10L) )
#'
#' # Show all parameters with defined values
#' pse$showParameters(  )
#'
#' # Show all parameters, including one that has an undefined value
#' # but an assigned deffinition (description)
#' pse$showParameters( na.rm=FALSE )
#' 
NULL

#' Self Variable Name
#'
#' Internal ParamSetI object method to determine what the object's
#' variable name is
#'
#' @name .selfVarName
#' @method .selfVarName ParamSetI
#' 
#' @details
#' 
#' \preformatted{
#' ## Method Usage:
#' myObject$.selfVarName( help=TRUE )
#'
#' myObject$.selfVarName( na.rm=TRUE )
#' }
#'
#' This is an internal method, you likely don't need to call it
#' yourself. It's designed to discover what the object variable name
#' is. This is then used in several places to produce non-generic
#' "real" examples for the user; The rationale is that they can then
#' copy-and-paste directly, without modifying a generic variable
#' name. Here it is used in \link{showParameters}.
#'
#' This turns out to be a REAL pain to do. I spent a lot of time
#' looking at deparse / sys.call to figure out what the "real"
#' variable name was. In particular this approach fails inside
#' implicit calls to show(). If show is explicitly called -
#' "show(myObject)" - then the variable name can be found on the call
#' stack. For some reason if myObject was printed (show()'n) implicitly
#' by calling just "myObject", then the stack
#' "does not go back far enough"
#'
#' So instead I am going to brute force look at the variable
#' space. I am only looking in the global environment (1), and
#' there are many ways this might not be quite right, but it
#' is close enough. I hope.
#' 
#' @param def Default 'myObj', the string to use if the method fails
#'     (finding the real object name is not always possible!)
#' @param fallbackVar Default "". Another default object name. To be
#'     honest, I forget why this was needed, but it
#'     was. For... reasons. I don't want to take it out to find out
#'     what breaks... It is used instead of \code{def}, when needed.
#' @param help Default FALSE. If TRUE, show this help and perform no
#'     other actions.
#'
#' @seealso \link{showParameters}
#' 
#' @examples
#'
#' ## This demo defines a toy object inheriting ParamSetI:
#' demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)
#' littleEngineThatCould <- ParamSetExample( params=list(inc=10L) )
#' littleEngineThatCould$.selfVarName()
#' 
NULL
