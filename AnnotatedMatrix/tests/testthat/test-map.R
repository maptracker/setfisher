library("AnnotatedMatrix")

symFile <- AnnotatedMatrix:::.makeTempFile("Symbol-To-Gene.mtx")
s2e <- AnnotatedMatrix(symFile)

test_that("Basic mapping", {
    qry <- c("GARS","SIGLEC19P")
    v1 <- s2e$map(qry)
    expect_identical(v1$Output,
                     c("LOC2617", "LOC2618", "LOC27036", "LOC114191"),
                     "Mapped output")
    expect_identical(v1$Score, c(4,2,2,1),
                     "Mapped score")
    expect_identical(v1$Factor,
                     c("Official", "Unofficial", "Unofficial", "Unknown"),
                     "Mapped factor levels")
    expect_identical(v1$Symbol, c("GARS", "GART", "SIGLEC7", "SIGLECP2"),
                     "Mapped metadata")
    expect_identical(attr(v1, "Mult.In"), c("GARS", "SIGLEC19P"),
                     "Mutliplicity attributes")


    v2 <- s2e$map(qry, format="vector")
    expect_equal(v2, c("LOC2617", "LOC2618", "LOC27036", "LOC114191"),
                 "Vector return", check.attributes = FALSE)
    expect_identical(names(v2), c("GARS", "GARS", "SIGLEC19P", "SIGLEC19P"),
                     "Vector names")
    
    
})

test_that("Collapsed mapping", {
    qry1 <- c("GARS","SIGLEC19P")
    v3 <- s2e$map(qry1, collapse="in")
    expect_identical(v3$Output,
                     c("LOC2617,LOC2618", "LOC27036,LOC114191"),
                     "Collapsed output")
    expect_identical(v3$Score, c(5,6),
                     "Collapsed score")
    expect_identical(v3$Factor,
                     c("Official,Unofficial", "Unofficial,Unknown"),
                     "Collapsed factor levels")

    qry2 <- c("P29","PRX", "GARS", "PAIS")
    v4 <- s2e$map(qry2, add.metadata=F, collapse='out')
    expect_identical(v4$Input,
                     c("P29", "P29,PRX", "P29", "PRX", "GARS",
                       "GARS,PAIS", "PAIS"),
                     "Collapsed input")
    expect_identical(v4$Output,
                     c("LOC5657", "LOC9588", "LOC25949", "LOC57716",
                       "LOC2617", "LOC2618", "LOC10606"),
                     "Collapsed output via input")
    
    
})
