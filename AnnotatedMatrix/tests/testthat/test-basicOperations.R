library("AnnotatedMatrix")

symFile <- AnnotatedMatrix:::.makeTempFile("Symbol-To-Gene.mtx")
test_that("Temp file copy", {
    expect_true(file.size(symFile) > 0,
                ".makeTempFile found sample file")    
})

s2e <- AnnotatedMatrix(symFile)
test_that("Object properties", {
    expect_identical("AnnotatedMatrix", class(s2e)[1],
                     "Proper class generated")
})

rn <- s2e$rNames()
cn <- s2e$cNames()

test_that("LoadingMTX", {
    expect_identical(167L, length(rn),
                     "Valid number of rows")
    expect_identical(67L, length(cn),
                     "Valid number of cols")
    expect_identical(208L, s2e$nnZero(),
                     "Valid number of connections")
    expect_identical(c("AR","AIRE1","IGF2RAS","IGF2R-AS1"), rn[c(1,2,166,167)],
                     "Spot check of rownames")
    expect_identical(c("LOC231","LOC326","LOC100188846","LOC100271873"),
                     cn[c(1,2,66,67)],
                     "Spot check of colnames")
    
})

test_that("Parameters", {
    expect_identical(s2e$param("name"), "Sample Human Symbol-to-Entrez Lookup" )
    expect_identical(s2e$param("rowdim"), "Gene Symbol"  )
    expect_identical(s2e$param("coldim"), "Entrez ID"  )
    expect_identical(s2e$param("cOlDiM"), "Entrez ID",
                     "case insensitivity")

    expect_identical(s2e$param("lemming"), NA,
                     "unset parameter")
    expect_identical(s2e$param("lemming", "fuzzy"), "fuzzy",
                     "setting parameter")
    expect_identical(s2e$param("lemming"), "fuzzy",
                     "setting parameter, confirming")

    val <- c("fuzzy", "wuzzy")
    expect_identical(s2e$param("lemming", val), "fuzzy",
                     "Default no vectors")
    
    expect_identical(s2e$param("lemming", val, is.scalar=FALSE), val,
                     "Setting a vector")
    expect_identical(s2e$param("lemming", "hippo", clobber=FALSE), val,
                     "No clobber")
})

test_that("Matrix recovery and transposition", {
    x <- s2e$matObj(transpose=TRUE)
    expect_identical(colnames(x), rn,
                     "Columns properly transposed")
    expect_identical(rownames(x), cn,
                     "Rows properly transposed")
    s2e$reset()
    nnz <- s2e$nnZero()
    s2e$filterByFactorLevel("Official")
    expect_true(s2e$nnZero() < nnz, "Connections drop after filtering")
    expect_true(s2e$nnZero(raw=TRUE) == nnz, "Raw object accessible")
    
})

test_that("Vector extraction from columns and rows", {
    s2e$reset()
    expect_identical(s2e$getRow("AIS"), c("LOC367", "LOC8626", "LOC260402"),
                     "Recover columns for row")
    expect_identical(s2e$getCol("LOC7038"), c("TG", "TGN", "AITD3"),
                     "Recover columns for row")
    
})

test_that("matrix extraction from columns and rows", {
    s2e$reset()
    expect_identical(s2e$getRow("AIS"), c("LOC367", "LOC8626", "LOC260402"),
                     "Recover columns for row")
    expect_identical(s2e$getCol("LOC7038"), c("TG", "TGN", "AITD3"),
                     "Recover columns for row")
    
    mChk <- matrix(c(2,2), ncol=1,
                   dimnames=(list("Entrez ID" = c("LOC7128","LOC28935"),
                                  "Gene Symbol" = c("A20"))))
    expect_identical(mChk, as.matrix(s2e$getRow("A20", format='matrix')),
                     "Matrix recovery of a row")
})

