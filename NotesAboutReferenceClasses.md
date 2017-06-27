## Adventures with R's ReferenceClass

#### Roxygen

The main shortcoming I've found is the ability to document parameters
for object methods. This is a [requested feature][Roxygen-196] that is
currently (Jun 2017) [being discussed][Roxygen-388] but not yet
implemented. I am currently putting in a docstring with a preformat
block to provide details on each method, eg:

```R
ParamSetI$methods(
    param = function (key = NA, val = NA, append = FALSE) {
        "\\preformatted{Gets/Sets a parameter. Parameters:
      key - The name of the parameter. If NA will simply return NA. Key names
            are case insensitive.
      val - Optional new value. If not NA, then the parameter value will be
            set, with behavior modified by some of the flags below
   append - Default FALSE, which will cause the parameter to be set to 'val'.
            If TRUE, val will be appended to the end of the current vector
            holding the value for 'key'.
}"
        if (is.na(key)) return( NA )
# ... etc more code ...

```

Within R this generates more whitespace (blank lines) than I'd like,
but it is mostly readable, both in the documentation and in the code.

#### Inheritance

`TL;DR`: Use the fully specified function name, eg `foo::bar` to
access `bar()` in package `foo`. Needed because even if Class X is
importing `foo`, when Class Y inherits (contains) X, foo is not included.

__Details:__ Imports do not appear to natively follow inherited
ReferenceClass classes. That is, say you have a class that is
importing `crayon` - a text coloring package that will be bringing in
the function `red()` - in this example:

```R
#' My Utilities
#' Does handy things
#' @import crayon
#' @importFrom methods new setRefClass
#' @export myUtils
#' @exportClass myUtils

myUtils <- setRefClass("foo",
  fields   = list( statVal = character ))
  
myUtils$methods(
  initialize = function(status="All clear", ...) {
     callSuper(...)
     statVal <<- status
  },
  showStatus = function() { message("Status is: ",red(statVal)) }
  )
```

This class will work fine on its own:

```R
> library("myUtils")
> mu <- myUtils()
> mu$showStatus()
Status is: All clear
```

If, however, you make another class that inherits the above one:

```R
#' My Tool
#' Does useful things
#' @import
#' @importFrom methods new setRefClass
#' @export myTool
#' @exportClass myTool

myTool <- setRefClass("foo",
  fields   = list( x = 'numeric' ),
  contains = c("myUtils"))
  
myTool$methods(
  initialize = function(status="SystemReady", x=42, ...) {
     callSuper(status=status, ...)
     x <<- x
  },
  showX = function() { showStatus(); message("x = ",x) }
  )
```

... then - presuming these classes are compiled and loaded as formal
packages and you do not have `crayon` loaded in your search space -
attempts to use myTool will throw an error:

```R
> library("myTool")
> mt <- myTool(x=99)
> mt$showX()
Error in showX() : could not find function "red"
```

I have spent a rather large amount of time trying to resolve this
problem. I have been able to address it by using `crayon::red()`
rather than `red()` in the code. I believe this should be a relatively
safe mechanism to handle the problem, but I am suspicious there are
more elegant approaches that I am not aware of.

This outcome apparently has to do with how ReferenceClasses handle
their environment; There is a section title _Inter-Package
Superclasses and External Methods_ in the `?ReferenceClasses`
documentation that might be referrencing the issue; However I had
trouble comprehending that section. It is possible

[Roxygen-196]: https://github.com/klutometis/roxygen/issues/196
[Roxygen-388]: https://github.com/klutometis/roxygen/issues/388
