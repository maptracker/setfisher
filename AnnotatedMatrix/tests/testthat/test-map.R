library("AnnotatedMatrix")

s2e <- AnnotatedMatrix( annotatedMatrixExampleFile() )

test_that("Basic mapping", {
    qry <- c("GARS","SIGLEC19P")
    v1 <- s2e$map(qry)
    expect_identical(v1$Output,
                     c("LOC2617", "LOC2618", "LOC27036"),
                     "Mapped output")
    expect_identical(v1$Score, c(4L,2L,2L),
                     "Mapped score")
    expect_identical(v1$Factor,
                     c("Official", "Unofficial", "Unofficial"),
                     "Mapped factor levels")
    expect_identical(v1$Symbol, c("GARS", "GART", "SIGLEC7"),
                     "Mapped metadata")
    expect_identical(rownames(v1),
                     c("GARS#1", "GARS#2", "SIGLEC19P"),
                     "Customized rownaming for non-unique entries")

    expect_identical(attr(v1, "Mult.In"), c("GARS"),
                     "Mutliplicity attributes")

    v2 <- s2e$map(qry, format="vector")
    expect_equal(unclass(v2), c("LOC2617", "LOC2618", "LOC27036"),
                 "Vector return", check.attributes = FALSE)
    expect_identical(names(v2), c("GARS", "GARS", "SIGLEC19P"),
                     "Vector names")
    
    
})

test_that("Collapsed mapping", {
    qry1 <- c("GARS","SIGLEC19P")
    v3 <- s2e$map(qry1, collapse="in")
    expect_identical(v3$Output,
                     c("LOC2617,LOC2618", "LOC27036"),
                     "Collapsed output")
    expect_identical(v3$Score, c(5L,2L),
                     "Collapsed score")
    expect_identical(v3$Factor,
                     c("Official,Unofficial", "Unofficial"),
                     "Collapsed factor levels")

    qry2 <- c("P29","PRX", "GARS", "PAIS")
    v4 <- s2e$map(qry2, add.metadata=FALSE, collapse='out')
    expect_identical(v4$Input,
                     c("P29", "P29,PRX", "P29", "PRX", "GARS",
                       "GARS,PAIS", "PAIS"),
                     "Collapsed input")
    expect_identical(v4$Output,
                     c("LOC5657", "LOC9588", "LOC25949", "LOC57716",
                       "LOC2617", "LOC2618", "LOC10606"),
                     "Collapsed output via input")
    
    
})

test_that("Appended Mapping", {
    myData <- data.frame(Flavor = c("Cherry", "Moose"),
                         Symbol = c("GARS","SIGLEC19P"),
                         stringsAsFactors=FALSE)
    v4 <- s2e$map( append.to=myData, append.col="Symbol", collapse="in" )

    addCol <- c("Output", "Score", "Factor", "Symbol.1", "Description")
    expect_identical(attr(v4, "Appended"), addCol,
                     "Check appended colnames")
    expect_identical(colnames(v4),
                     c("Flavor", "Symbol", addCol),
                     "Appended data frame")
    expect_identical(v4$Output, c("LOC2617,LOC2618", "LOC27036"),
                     "Collapsed output")
    
})
