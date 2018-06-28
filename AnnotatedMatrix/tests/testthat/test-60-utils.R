library("AnnotatedMatrix")

message("Testing: Utility functions")

test_that("Parameter management", {

    param <- "fooTest"
    lcp   <- tolower(param)
    def   <- "A default value"
    vOpt  <- "A value set by option"
    val   <- "Some value I set"
    expect_null(annotatedMatrixParameter(param),
                info="Un-set parameter should be NULL initially")

    namedDef <- stats::setNames(def, param)
    expect_identical(annotatedMatrixParameter(param, default=def),
                     namedDef, info="Define a default value")


    ## Set an option with the key defined by a variable:
    optArgs <- list()
    ## First, do it without the "annotatedmatrix" prefix:
    optArgs[[ lcp ]] <- vOpt
    do.call("options", optArgs)
    expect_null(annotatedMatrixParameter(param),
                info="Still should be unfound")
    ## Now set with the proper prefix:
    optArgs[[ paste0("annotatedmatrix", lcp) ]] <- vOpt
    do.call("options", optArgs)
    namedOpt <- stats::setNames(vOpt, param)
    expect_identical(annotatedMatrixParameter(param), namedOpt,
                     info="Value recovered from option fallback")


    ## Now set it explicitly:
    namedVal <- stats::setNames(val, param)
    expect_identical(annotatedMatrixParameter(param, val), namedVal,
                     info="Value recovered from explicit set")
    ## Should "still be there" afterwards:
    expect_identical(annotatedMatrixParameter(param), namedVal,
                     info="Value persists after explicit setting")

    ## Non-vector values
    myDF <- data.frame(x =1:3, y=4:6)
    dfP  <- "CoolDataFrame"
    expect_identical(annotatedMatrixParameter(dfP, myDF), myDF,
                     info="Non-vector parameters should be returned un-named")

    ## Multiple values
    x <- annotatedMatrixParameter(c(dfP, param))
    expect_identical(storage.mode(x), "list",
                     info="Multiple values returned as list by default")

    expect_identical(x[[ param ]], val,
                     info="Check list names [[1]]")
    expect_identical(x[[ dfP ]], myDF,
                     info="Check list names [[2]]")
    
    
 })
