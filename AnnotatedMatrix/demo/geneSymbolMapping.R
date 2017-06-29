## Work from a copy of a sample file provided in inst/
library("AnnotatedMatrix")
src <- file.path(path.package("AnnotatedMatrix"),
                 "extdata/Symbol-To-Gene.mtx")
## (for developing - files will still be in inst/ subdirectory)
if (!file.exists(src)) src <- file.path(path.package("AnnotatedMatrix"),
                                        "inst/extdata/Symbol-To-Gene.mtx")
trg <- file.path(tempdir(), "Symbol-To-Gene.mtx")
file.copy(src,trg)
message("Working with MatrixMarket file:\n    file:\n  ",
        trg, "\n    from: ",src)

s2e <- AnnotatedMatrix(trg)

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
