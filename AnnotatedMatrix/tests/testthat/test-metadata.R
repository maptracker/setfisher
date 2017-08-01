library("AnnotatedMatrix")

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )

test_that("Basic mapping", {
    id  <- c("LOC8784", "LOC25824")
    key <- c("Symbol")

    ## No query, all metadata
    allMD <- s2e$metadata()
    expect_identical(length(rownames(allMD)), 67L,
                     "Full metadata recovery, row count")
    expect_identical(length(colnames(allMD)), 3L,
                     "Full metadata recovery, col count")

    ## Query by both ID and key
    both <- s2e$metadata(id, key)
    bothExpct <- setNames(c("TNFRSF18", "PRDX5"), id)
    expect_identical(both, bothExpct,
                     "Named vector for single column")
    
    ## Query by ID
    idOnly <- s2e$metadata(id)
    expect_identical(idOnly$Symbol, c("TNFRSF18", "PRDX5"),
                     "ID-only recovery")
    expect_identical(idOnly$id, id,
                     "ID-only recovery (query check)")

    ## Query by Key
    keyOnlyVec <- s2e$metadata(key=key)
    ## Row order can change in data.table objects. Sort the result 
    keyOnlyVec <- keyOnlyVec[ order(keyOnlyVec) ]

    expect_true(is.vector(keyOnlyVec), "Single column, vector recovery")
    expect_identical(length(keyOnlyVec), 67L, "Single column, vector recovery")
    expect_equivalent(head( keyOnlyVec, 4 ), c("AIPL1", "AIR", "AIRE", "AIRN"),
                      "Key-only, single column, value check")
    expect_identical(head( names(keyOnlyVec), 4 ),
                     c("LOC23746", "LOC7808", "LOC326", "LOC100271873"),
                     "Key-only, single column, name check")
})

test_that("Unknown queries", {
    id  <- c("LOC8784", "LOC25824")
    key <- c("Symbol")
    
    ## Query by both ID and key plus bogus key
    expect_message(bogus <- s2e$metadata(id, c("Symbol", "Pirate")))
    bothExpct <- setNames(c("TNFRSF18", "PRDX5"), id)
    expect_identical(bogus, bothExpct,
                     "Named vector for single column")
    
    ## Empty columns
    silly <- c("Humpty", "Dumpty")

    ## Columns kept if both id and key are set:
    notThere <- s2e$metadata(id=silly, key=key)
    ntRv <- setNames(as.character(c(NA,NA)), silly)
    expect_identical(notThere, ntRv)
    stillNotThere <- s2e$metadata(id=silly, key=key, na.rm=TRUE)
    expect_identical(stillNotThere, ntRv)
    
      
})
