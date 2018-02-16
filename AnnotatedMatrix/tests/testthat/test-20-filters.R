library("AnnotatedMatrix")

message("Testing: Filters")

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )

test_that("Matrix filters", {
    afCon <- 180L
    expect_identical(afCon, s2e$nnZero(),
                     info="Valid number of auto-filtered connections")
    
    s2e$reset()
    rawCon <- 208L
    expect_identical(rawCon, s2e$nnZero(),
                     info="Valid number of connections after reset")

    ## Score filters
    s2e$filterByScore(min=2, max=3)
    expect_identical(152L, s2e$nnZero(),
                     info="Valid number of connections after score filter")
    s2e$reset()
    expect_null(s2e$matrixUse, "Clearing of $matrixUse after reset")
    expect_identical(rawCon, s2e$nnZero(),
                     info="Valid number of connections after reset")
    
    ## Filter by factor levels
    s2e$reset()
    expect_identical(s2e$map("LOC7038")$Output, c("TG", "TGN", "AITD3"),
                     info="Unfiltered data")
    s2e$filterByFactorLevel("Official")
    
    expect_identical(s2e$map("LOC7038")$Output, c("TG"),
                     info="Official symbols only")

    ## Filter by metadata text
    s2e$reset()
    sym <- "AITD3"
    expect_equivalent(unclass(s2e$map(sym, format="vector")),
                      c("LOC7038", "LOC57623", "LOC387580"),
                      info="Symbol with scruffy link")
    ## Remove genes with "{Deprecated}" in the description:
    s2e$filterByMetadata(key="Description", keep=FALSE, val="{Deprecated}")
    expect_equivalent(unclass(s2e$map(sym, format="vector")),
                      c("LOC7038", "LOC57623"),
                      info="Symbol with scruffy link, unscrufified")

    ## Check summary information
    fs <- s2e$filterSummary()
    expect_identical(fs[1, "metric"], "Description LIKE {Deprecated}",
                     info="Summary metric")
    expect_identical(fs[1, "Col"], 5L,
                     info="Summary counts")

    ## Filter by row and column names
    s2e$reset()
    bogus  <- "FlamingMonkeys"
    newRow <- c("ZFAT", bogus, "TGN")
    s2e$rNames( newRow )
    expect_identical(s2e$nnZero(), 2L, "Restricted rows")
    ## Some issues with Matrix falling back to dgCMatrix:
    expect_identical(class(s2e$matrixUse)[1], "dgTMatrix",
                     info="Class safety check")
    expect_identical(s2e$map(bogus)[bogus, "Score"], 0L,
                     info="Verify novel row is recognized (Score = 0 != NA)")

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
    
    ## Filter by ID
    s2e$reset()
    rN <- s2e$rNames()
    expect_identical(length(rN), 167L, "Sanity check")
    keepRow <- c("AR","AIRE1","APSI","PGA1")
    crc <- s2e$filterById(keepRow, keep=TRUE, MARGIN='row')
    expect_equivalent(crc, c(202L, 163L, 63L),
                      info="Heavy filter keeping 4 IDs")
    expect_identical(s2e$rNames(nonzero=TRUE), keepRow,
                     info="Filter keeps just 4 IDs")
    s2e$reset()
    ## Restrict to some weirdly specific symbol patterns.
    s2e$filterById(c('[a-z][0-9][a-z][0-9][a-z]','[0-9]-[0-9]'),
                   keep=TRUE, exact=FALSE, MARGIN='row')
    expect_identical(s2e$rNames(nonzero=TRUE),
                     c("RP3-426F10.1", "RP11-397D12.1", "NT5C1B"),
                     info="Regular expression ID filter")

    s2e$reset()
    ## Verify that we get a warning when we remove everything when
    ## MARGIN is NULL
    expect_message(s2e$filterById("AIRE1", keep=TRUE),
                   info="removed ALL matrix entries")
    
    ## Applied filters field
    s2e$reset()
    s2e$autoFilter()
    af <- s2e$appliedFilters()
    sf <- s2e$setFilters
    expect_identical(af, sf,
                     info="Method without args should just be setFilters")
    expect_identical(length(af), 2L, "Two automatic filters")
    expect_identical(af[1], "LEVELS == Unknown ## Unknown status = uncertain provenance in MapTracker", info="First example filter")

    
   
})
