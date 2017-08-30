library("AnnotatedMatrix")

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )

test_that("Matrix filters", {
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
    expect_identical(s2e$map("LOC7038")$Output, c("TG", "TGN", "AITD3"),
                     "Unfiltered data")
    s2e$filterByFactorLevel("Official")
    
    expect_identical(s2e$map("LOC7038")$Output, c("TG"),
                     "Official symbols only")

    ## Filter by metadata text
    s2e$reset()
    sym <- "AITD3"
    expect_equivalent(unclass(s2e$map(sym, format="vector")),
                      c("LOC7038", "LOC57623", "LOC387580"),
                      "Symbol with scruffy link")
    ## Remove genes with "{Deprecated}" in the description:
    s2e$filterByMetadata(key="Description", keep=FALSE, val="{Deprecated}")
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
    expect_identical(s2e$map(bogus)[bogus, "Score"], 0L,
                     "Verify novel row is recognized (Score = 0 != NA)")

    ## Filter by row or column counts
    s2e$reset()
    rc <- s2e$rCounts()
    expect_identical(sum( rc > 1 ), 24L, "Sanity check on row counts")
    rfilt <- s2e$filterByCount("row", min=2, filterEmpty=TRUE,
                               reason="Keep only symbols with 2+ genes")
    expect_identical(unname(rfilt), c(143L,143L,18L), "Row Count filters")
    expect_identical(nrow(s2e$matrixUse), 24L, "Row Count filters")
    cfilt <- s2e$filterByCount("col", min=3, max=5, filterEmpty=TRUE,
                               reason="Genes with 3-5 symbols")
    expect_identical(unname(cfilt), c(62L,21L,48L), "Col Count filters")
    ## Only a few entries left
    expect_identical(s2e$rNames(), c("p75", "SIGLEC19P", "SIGLECP2"))
    expect_identical(s2e$cNames(), c("LOC27036"))
    
    
    
    ## Applied filters field
    s2e$reset()
    s2e$autoFilter()
    af <- s2e$appliedFilters()
    sf <- s2e$setFilters
    expect_identical(af, sf, "Method without args should just be setFilters")
    expect_identical(length(af), 2L, "Two automatic filters")
    expect_identical(af[1], "LEVELS == Unknown ## Unknown status = uncertain provenance in MapTracker", "First example filter")

    
   
})
