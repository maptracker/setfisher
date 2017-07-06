library("AnnotatedMatrix")

## Load a small matrix that provides a lookup of gene symbols to their
## Entrez Gene accessions.
## ( .makeTempFile() is a utility method that's only of use in these demos.
##   Normally you would be providing your own static path to the file. )
symFile <- AnnotatedMatrix:::.makeTempFile("Symbol-To-Gene.mtx")
s2e     <- AnnotatedMatrix(symFile)

## -----------------------------------------------------------------------##
s2e$message("Handling multiplicity by concatenation", bgcolor='cyan')

## Here we are recovering symbols for two gene accessions:
s2e$map(c("LOC7038","LOC51248"))

## To keep the output with one row per gene, we can just collapse the
## symbols into a single token-separated string:
s2e$map(c("LOC7038","LOC51248"), collapse="in")

## We can configure how multiple cells are combined into one:
s2e$map(c("LOC7038","LOC51248"), collapse="in",
        collapse.token=' | ', # The text used to separate IDs
        collapse.factor=max ) # The function used to aggregate scores

## -----------------------------------------------------------------------##
s2e$message("Dealing with well-behaved maps", bgcolor='cyan')

## Sometimes you have mappings that are not unique, but are
## differentially scored to highlight the 'better' mappings:
s2e$map(c("AR","AIS1"))

## In these cases, the keep.best parameter can be used to only keep
## the best mapped output(s) for each input. Here, the two symbols
## each have a unique Official symbol, and our results can be
## collapsed to unique rows with keep.best
s2e$map(c("AR","AIS1"), collapse='in', keep.best=TRUE)

## -----------------------------------------------------------------------##
s2e$message("When you absolutely must get to just one row", bgcolor='cyan')

## The symbol "p40" is not official, and it is an unofficial symbol
## for a LOT of loci:
s2e$map("p40")

## From the matrix, we really do not have a good mechanism to pick one
## of these mappings over the others. If we really MUST have a unique
## mapping, we can use the provided takeLowestThing() function to pick
## the Entrez accession with the smallest integer:
s2e$map("p40", collapse="input", collapse.name=takeLowestThing)

####################################################################
## Note that this is an extremely crude way to resolve the issue! ##
## Ideally, attempt to recover the data with better resolved      ##
## identifiers. Or get a 1:1 mapping from the original source.    ##
## Ultimately, use such hacks only when aboslutely neccesary and  ##
## when all other recourses have failed.                          ##
####################################################################
