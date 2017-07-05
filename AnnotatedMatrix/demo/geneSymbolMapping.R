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

symFile <- tmpFile("Symbol-To-Gene.mtx")
s2e <- AnnotatedMatrix(symFile)

mySyms <- unlist(strsplit("AIRE GITR-D p40 p75 SHREW-1 NSGPx BorkBork AOP2 PRX", ' '))

## Map those symbols over to genes
s2e$actionMessage("Mapping symbols to genes")
genes1 <- s2e$map(mySyms)
print(genes1)

## Messy! Try keeping the best 
s2e$actionMessage("Mapping symbols to genes, keeping only best score")
genes2 <- s2e$map(mySyms, keep.best=TRUE)
print(genes2)

## Did not really help! Collapse by input ID:
s2e$actionMessage("Mapping symbols to genes, one row per input")
genes3 <- s2e$map(mySyms, collapse='in', keep.best=TRUE)
print(genes3)

## Still ugly! The Unofficial symbols are irritating, let's remove them
s2e$actionMessage("Mapping symbols to genes, only official and preferred symbols")
s2e$filterByScore(min=3)
genes4 <- s2e$map(mySyms, collapse='in', keep.best=TRUE)
print(genes4)

## Bad idea, most of those symbols are unofficial!
## Reset the filters:
s2e$reset()
## For symbols with multiple mappings, use the built-in desperation
## method that picks "the one with the lowest number"
genes5 <- s2e$map(mySyms, collapse='in', keep.best=TRUE,
                  collapse.name=takeLowestThing)
print(genes5)
