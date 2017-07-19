library("ParamSetI")
## Set up a toy class to inherit the module:

foo <- setRefClass("foo",
  fields   = list( x = 'integer' ),
  contains = c("ParamSetI"))


foo$methods(
    initialize = function( x=1L, params=NULL, ... ) {
        x <<- x
        ## Set some default parameters, one for a key with a
        ## definition, the other without one
        setParamList(list(myNum=2.4, bar=3L))
        ## Then set any user values based on params
        callSuper(params=params, ...)
        defineParameters("
myNum  [numeric] A number for me
myInt  [integer] An integer for me
myChar [character] A string for me
")
    })

psObj <- foo()
psObj2 <- foo( params=list(color='navy', size='medium') )

test_that("Object creation", {
    
    expect_identical(psObj$allParams(), c("mynum", "myint", "mychar", "bar"),
                     "Check recognized params")
    
    expect_identical(psObj$param("MYNUM"), 2.4, "case-insensitive")
    expect_identical(psObj$param("bar"), 3L, "default, not defined")
    expect_identical(psObj$param("lemming"), NA, "unset")
    
    expect_true(psObj$hasParam("myNum"), "Described and set")
    expect_true(psObj$hasParam("myInt"), "Described, not set")
    expect_false(psObj$hasParam("foodles"), "Neither set nor described")

    ## Definitions and Class restrictions
    expect_equivalent(psObj$paramDefinition('mynum'), "A number for me",
                      "Definition recovery")
    expect_equivalent(psObj$paramClass('mynum'), "numeric",
                      "Definition recovery")

    ## Name manipulation
    expect_equivalent(psObj$paramName('mynum'), "myNum",
                      "Case-preserved name recovery")
    expect_equivalent(psObj$paramName('mynum', "MyNum"), "MyNum",
                      "Change case of name")
    expect_equivalent(psObj$paramName('mynum', "CarpeDiem"), "MyNum",
                      "Can't change the name of, uh, name")

    ## Passing values on object creation
    expect_identical(psObj2$hasParam(c("color", "size")), c(TRUE,TRUE),
                      "New params on init")
    expect_identical(psObj2$param("size"), "medium",
                      "New params on init")
    
    
})
                   
test_that("Parameter setting", {
    oldKey1 <- "myNum"
    newVal1 <- 12.34
    expect_identical(psObj$param(oldKey1, newVal1), newVal1, "Reset value")
    expect_identical(psObj$param(oldKey1), newVal1, "new value persistent")
    
    expect_identical(psObj$param(oldKey1, 9999, clobber=FALSE), newVal1,
                     "clobber protects existing values")

    newKey1 <- "flower"
    newVal2 <- "pansy"
    expect_identical(psObj$param(newKey1, newVal2, clobber=FALSE), newVal2,
                     "clobber irrelevant for new keys")

    
    

    vecVal <- c(1.2, 2.3, 3.4)
    newKey2 <- "someScalar"
    psObj$param(newKey2, vecVal)
    expect_identical(psObj$param(newKey2), vecVal[1], "default to scalar")
    addVal <- 5.6
    psObj$param(newKey2, addVal, append=TRUE)
    expect_identical(psObj$param(newKey2), c(vecVal[1], addVal), "appending")
    
    newKey3 <- "someVector"
    psObj$param(newKey3, vecVal, is.scalar=FALSE)
    expect_identical(psObj$param(newKey3), vecVal, "is.scalar override")
    
    
    oldKey2 <- "myInt"
    notInt  <- 3.14
    asInt   <- as.integer(notInt)
    expect_identical( psObj$param(oldKey2, notInt), asInt, "Default coercion")
    expect_identical( psObj$param(oldKey2, notInt, coerce=FALSE), NA,
                     "Coercion turned off")
    expect_identical( psObj$param(oldKey2, notInt, coerce=FALSE,
                                  check.class=FALSE), notInt,
                     "Coercion and class checking turned off")

    ## Set by list
    dVec <- c("bears","beets","BSG")
    psObj$setParamList(list(x=17, y=42, dwight=dVec))
    expect_identical( psObj$param("x"), 17,
                     "Set values by list 1")
    expect_true( psObj$hasParam("y") )
    expect_identical( psObj$param("dwight"), dVec[1],
                     "Set values by list 2 (forced scalar)")
    psObj$setParamList(list(x=17, y=42, dwight=dVec), is.scalar=FALSE)
    expect_identical( psObj$param("dwight"), dVec,
                     "Set values by list 2 (vector)")
})


test_that("EventLogger Inheritance", {
    expect_message( psObj$message("This is a message") )
    expect_message( psObj$err("This is an error message") )
    expect_message( psObj$dateMessage("This is a dated message") )
})
