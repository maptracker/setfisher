library("AnnotatedMatrix")

symFile <- AnnotatedMatrix:::.makeTempFile("Symbol-To-Gene.mtx")
s2e <- AnnotatedMatrix(symFile)

test_that("LoadingMTX", {
    rawCon <- 208L
    s2e$reset()
    expect_identical(rawCon, s2e$nnZero(),
                     "Valid number of connections")
    s2e$filterByScore(min=2, max=3)
    expect_identical(152L, s2e$nnZero(),
                     "Valid number of connections after score filter")
    s2e$reset()
    expect_null(s2e$matrixUse, "Clearing of $matrixUse after reset")
    expect_identical(rawCon, s2e$nnZero(),
                     "Valid number of connections after reset")
    
    
    expect_identical(s2e$getCol("LOC7038"), c("TG", "TGN", "AITD3"),
                     "Unfiltered data")
    s2e$filterByFactorLevel("Official")
    
    expect_identical(s2e$getCol("LOC7038"), c("TG"),
                     "Official symbols only")

                     
    
})
