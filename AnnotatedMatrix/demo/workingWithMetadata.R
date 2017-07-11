library("AnnotatedMatrix")

## .makeTempFile() is a utility for handling demo files - you normally
## would pass the path to your matrix file directly

lolFile <- AnnotatedMatrix:::.makeTempFile("ListOfLists.inp")
lol     <- AnnotatedMatrix(lolFile)

lol$message("Loaded List-of-list sample file", bgcolor='cyan')
lol$actionMessage("The file contains embeded metadata, plus two sidecars")


lol$actionMessage("Get the codenames for all campers:")
lol$metadata(key="CodeName")

lol$actionMessage("Find the inexpensive campers:")
rate <- lol$metadata(key='Rate')
rate[ rate < 30 ]
