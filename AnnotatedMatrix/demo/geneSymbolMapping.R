## Work from a copy of a sample file provided in inst/
library("AnnotatedMatrix")
src <- file.path(path.package("AnnotatedMatrix"),
                 "extdata/Symbol-To-Gene.mtx")
## (for developing - files will still be in inst/ subdirectory)
if (!file.exists(src)) src <- file.path(path.package("AnnotatedMatrix"),
                                        "inst/extdata/Symbol-To-Gene.mtx")
trg <- file.path(tempdir(), "Symbol-To-Gene.mtx")
file.copy(src,trg)
message("Working with MatrixMarket file: ", trg, " from ",src)

s2e <- AnnotatedMatrix(trg)
