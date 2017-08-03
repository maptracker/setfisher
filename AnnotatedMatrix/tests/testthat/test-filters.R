library("AnnotatedMatrix")

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )

test_that("LoadingMTX", {
    afCon <- 180L
    expect_identical(afCon, s2e$nnZero(),
                     "Valid number of auto-filtered connections")
    
    s2e$reset()
    rawCon <- 208L
    expect_identical(rawCon, s2e$nnZero(),
                     "Valid number of connections after reset")

    ## Score filters
    s2e$filterByScore(min=2, max=3)
    expect_identical(152L, s2e$nnZero(),
                     "Valid number of connections after score filter")
    s2e$reset()
    expect_null(s2e$matrixUse, "Clearing of $matrixUse after reset")
    expect_identical(rawCon, s2e$nnZero(),
                     "Valid number of connections after reset")
    
    ## Filter by factor levels
    s2e$reset()
    expect_identical(s2e$getCol("LOC7038"), c("TG", "TGN", "AITD3"),
                     "Unfiltered data")
    s2e$filterByFactorLevel("Official")
    
    expect_identical(s2e$getCol("LOC7038"), c("TG"),
                     "Official symbols only")

    ## Filter by metadata text
    s2e$reset()
    sym <- "AITD3"
    expect_equivalent(unclass(s2e$map(sym, format="vector")),
                      c("LOC7038", "LOC57623", "LOC387580"),
                      "Symbol with scruffy link")
    ## Remove genes with "{Deprecated}" in the description:
    s2e$filterByMetadata(key="Description", val="{Deprecated}")
    expect_equivalent(unclass(s2e$map(sym, format="vector")),
                      c("LOC7038", "LOC57623"),
                      "Symbol with scruffy link, unscrufified")

    ## Check summary information
    fs <- s2e$filterSummary()
    expect_identical(fs[1, "metric"], "Description LIKE {Deprecated}",
                     "Summary metric")
    expect_identical(fs[1, "Col"], 5L,
                     "Summary counts")

    ## Filter by row and column names
    s2e$reset()
    bogus  <- "FlamingMonkeys"
    newRow <- c("ZFAT", bogus, "TGN")
    s2e$rNames( newRow )
    expect_identical(s2e$nnZero(), 2L, "Restricted rows")
    ## Some issues with Matrix falling back to dgCMatrix:
    expect_identical(class(s2e$matrixUse)[1], "dgTMatrix", "Class safety check")
    expect_identical(s2e$map(bogus)[bogus, "Score"], 0,
                     "Verify novel row is recognized (Score = 0 != NA)")
    
})
