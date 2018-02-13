library("AnnotatedMatrix")

message("Testing: $product()")

test_that("Product of two matrices", {
    ## A Pair of toy matrices with overlap
    mach <- AnnotatedMatrix( annotatedMatrixExampleFile("Machines.gmt") )
    comp <- AnnotatedMatrix( annotatedMatrixExampleFile("Components.gmt") )
    expect_message(prod <- mach$product(comp),
                   "You should really specify this value yourself",
                   "Warning if valfunc is not set")

    ## Confirm that the transitive mappings worked
    tmap <- list(Collector=c("Bolt", "Gear", "Rod", "SteelPlate"),
                 RockDrill=c("Bolt", "Gear", "Rod", "SteelPlate", "Blade"),
                 Rover=c("Gear", "Rod", "SteelPlate"))
    for (m in names(tmap)) {
        got  <- sort(prod$map(m)$Output)
        want <- sort(tmap[[ m ]])
        expect_identical(got, want, paste("Expected transitive map for", m))
    }

    ## Check that metadata carried over:
    md <- list('Actuator'="Produces linear motion",
               'Bumper'="Collision protection",
               'Console'="Example with no connections to Machines.gmt",
               'Crystal'=as.character(NA),
               'CuttingBit'="For boring holes",
               'DriveShaft'="Transfers Power",
               'GearBox'="For changing speed",
               'MountingPlate'="For attaching to surfaces",
               'Wormhole'=as.character(NA))
    mdKey <- "Description"
    for (x in names(md)) {
        got  <- prod$metadata(x, mdKey)
        want <- md[[x]]
        expect_equivalent(got, want, paste("Expected", mdKey, "for", x))
    }

})
