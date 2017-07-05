## Work from a copy of a sample file provided in inst/
library("AnnotatedMatrix")
srcDir <- path.package("AnnotatedMatrix")
## (for developing - files will still be in inst/ subdirectory):
if (file.exists(file.path(srcDir, "inst"))) srcDir <- file.path(srcDir, "inst")
srcDir <- file.path(srcDir, "extdata")
tmpDir <- tempdir()
tmpFile <- function(name) {
    src    <- file.path(srcDir, name)
    trg    <- file.path(tmpDir, name)
    file.copy(src,trg)
    sfx    <- gsub('.+\\.', '', name)
    message("Working with ", sfx, " file:\n      copy: ",
            trg, "\n    source: ",src)
    trg
}

rcpFile <- tmpFile("Recipes.gmt")
rcp     <- AnnotatedMatrix(rcpFile)

rcp$message("Recipes using potato:", color='red')
print(rcp$map("Potato")$Output)

rcp$message("Potato recipes, DF collapsed to unique inputs:", color='red')
print(rcp$map("Potato", collapse='Input'))

rcp$message("Ingredients in Baklava:", color='red')
print(rcp$map("baklava")$Output)


