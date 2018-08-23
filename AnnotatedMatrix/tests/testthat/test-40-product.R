library("AnnotatedMatrix")

message("Testing: $product()")

test_that("Identifying shared dimensions", {
    ## A Pair of toy matrices with overlap
    mach <- AnnotatedMatrix( annotatedMatrixExampleFile("Machines.gmt") )
    comp <- AnnotatedMatrix( annotatedMatrixExampleFile("Components.gmt") )

    sd    <- mach$sharedDimensions(comp)
    expSd <- stats::setNames(c(2L, 1L), c("Lft","Rgt"))
    
    expect_identical(sd, expSd, info="Basic shared dimension calculation")
    
    sd   <- mach$sharedDimensions(comp, dim2=1)
    expect_identical(sd, expSd, info="Manual specification with numeric 1")
    
    sd   <- mach$sharedDimensions(comp, dim1=2)
    expect_identical(sd, expSd, info="Manual specification with numeric 2")
    
    sd   <- mach$sharedDimensions(comp, dim2='1')
    expect_identical(sd, expSd, info="Manual specification with '1'")
    
    sd   <- mach$sharedDimensions(comp, dim1='2')
    expect_identical(sd, expSd, info="Manual specification with '2'")
    
    sd   <- mach$sharedDimensions(comp, dim2='row')
    expect_identical(sd, expSd, info="Manual specification with 'row'")
    
    sd   <- mach$sharedDimensions(comp, dim1='col')
    expect_identical(sd, expSd, info="Manual specification with 'col'")

    ## Check that we can silence verbosity:
    expect_silent(sd <- mach$sharedDimensions(comp, dim1='row', warn=FALSE))
    
    expect_identical(sd, stats::setNames(as.integer(c(NA, NA)), c("Lft","Rgt")),
                     info="Default failure is NA dimensions")

    expect_message(sd <- mach$sharedDimensions(comp, dim1='row', fail.na=FALSE),
                   "Matrices do not appear to have a shared dimension",
                   info="Verify that warning is shown by default")

    expect_true(sd[1] > 0,
                info="If fail.na is FALSE, should still get dimension")
    
})

test_that("Product of two matrices", {
    ## A Pair of toy matrices with overlap
    mach <- AnnotatedMatrix( annotatedMatrixExampleFile("Machines.gmt") )
    comp <- AnnotatedMatrix( annotatedMatrixExampleFile("Components.gmt") )
    expect_message(prod <- mach$product(comp),
                   "You should really specify this value yourself",
                   info="Warning if valfunc is not set")

    ## Confirm that the transitive mappings worked
    tmap <- list(Collector=c("Bolt", "Gear", "Rod", "SteelPlate"),
                 RockDrill=c("Bolt", "Gear", "Rod", "SteelPlate", "Blade"),
                 Rover=c("Gear", "Rod", "SteelPlate"))
    for (m in names(tmap)) {
        got  <- sort(prod$map(m)$Output)
        want <- sort(tmap[[ m ]])
        expect_identical(got, want,
                         info=paste("Expected transitive map for", m))
    }

    ## Use IJX forms of these toy matrices to check product score
    ## mapping more closely

    mach <- AnnotatedMatrix( annotatedMatrixExampleFile("Machines.ijx") )
    comp <- AnnotatedMatrix( annotatedMatrixExampleFile("Components.ijx") )

    ## `mach` tells us how many components are in a machine
    ## `comp` counts the number of parts that are in a component
    
    ## We want the transitive product that tells us how many parts are
    ## in a machine. In this case, this will be essentially the same
    ## as a 'traditional' matrix product - the sum of the products. We
    ## need to write a callback method for that:

    matProd <- function(l, r) sum( l * r )
    
    prod <- mach$product(comp, valfunc=matProd)
    ## Verify that metadata carried over properly:
    expect_identical(prod$rNames(), c("Collector", "RockDrill", "Rover"),
                     info="$product() rownames are correct")
    expect_identical(prod$cNames(), c("Bolt", "Gear", "Rod", "SteelPlate", "Blade"),
                     info="$product() colnames are correct")
    
    ## Inspect the actual matrix
    pMat <- prod$matObj()
    expect_identical(names(dimnames(pMat)), c("Machine", "Part"),
                     info="$product() dimensions properly transfered")

    ## Verify mapped counts
    partCount <- list(Collector=c(4, 8, 4, 1, 0),
                      RockDrill=c(24, 9, 1, 12, 1),
                      Rover=c(0, 7, 10, 18, 0))
    
    for (component in names(partCount)) {
        pc <- partCount[[ component ]]
        expect_identical(unname(pMat[component,]), pc,
                         info=paste("Expected part count for", component))
    }

    ## Check that metadata carried over:
    md <- list(Blade="Cutting part",
               Bolt="Threaded attachment",
               Gear="Rotational power transfer",
               Rod="Stiff cylinder",
               SteelPlate="A flat metal surface",
               
               Collector="Gathers samples",
               RockDrill="Extracts samples",
               Rover="Wheeled explorer")
    mdKey <- "Description"
    for (x in names(md)) {
        got  <- prod$metadata(x, mdKey)
        want <- md[[x]]
        expect_equivalent(got, want,
                          info=paste("Expected", mdKey, "for", x))
    }

### TODO - Check dimension names, metadata for the IJX loads
### TODO - Add parameter checking when it's implemented in AnnotatedMatrix.R
    
})
