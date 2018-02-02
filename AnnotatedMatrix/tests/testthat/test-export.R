library("AnnotatedMatrix")

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )

test_that("Export functions", {

    ## Filter the matrix so that it is smaller:
    s2e$rNames(c("HFH2","AIS1","p63"))

    mdf <- s2e$melt()
    ## Not sure we can count on the order of the rows, order them here
    mdf <- mdf[ order(mdf[[1]], mdf[[2]], mdf[[3]]), ]

    expect_identical(colnames(mdf), c("Symbol","EntrezGene", "Val"),
                     "Column name remapping")
    ## Verify DF contents:
    expect_identical(mdf[[1]], c("AIS1", "AIS1", "AIS1", "HFH2", "p63", "p63"))
    expect_identical(mdf[[2]], c("LOC260402", "LOC27022", "LOC373071", "LOC27022", "LOC10970", "LOC7405"))
    expect_identical(mdf[[3]], c(2, 2, 4, 2, 2, 2))

    gmt <- s2e$as.gmt()
    expect_identical(s2e$as.gmt(),
                     c("HFH2\tNA\tLOC27022\n",
                       "AIS1\tNA\tLOC27022\tLOC260402\tLOC373071\n",
                       "p63\tNA\tLOC7405\tLOC10970\n"), "GMT format")
   
})
