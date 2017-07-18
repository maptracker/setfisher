library("AnnotatedMatrix")

## Load a small silly GMT file that relates some foods to their ingredients
## ( .makeTempFile() is a utility method that's only of use in these demos.
##   Normally you would be providing your own static path to the file. )


### ---------------------------------------------------- ###
rcpFile <- AnnotatedMatrix:::.makeTempFile("Recipes.gmt")
rcp     <- AnnotatedMatrix(rcpFile)
rcp$message("Loaded GMT sample file", bgcolor='cyan')

rcp$message("Recipes using potato:", color='yellow')
print(rcp$map("Potato", warn=FALSE)$Output)

rcp$message("Potato recipes, DF collapsed to unique inputs:", color='yellow')
print(rcp$map("Potato", collapse='Input'))

rcp$message("Ingredients in Baklava:", color='yellow')
print(rcp$map("baklava", warn=FALSE)$Output)


### ---------------------------------------------------- ###
lolFile <- AnnotatedMatrix:::.makeTempFile("ListOfLists.inp")
lol     <- AnnotatedMatrix(lolFile)
lol$message("Loaded List-of-list sample file", bgcolor='cyan')

lol$message("People who are both chefs and managers:", color='yellow')
print(lol$map(c("Chefs", "Managers"), collapse='out'))

lol$message("Things Henry can do:", color='yellow')
print(lol$map("Henry"))


### ---------------------------------------------------- ###
dirPath <- AnnotatedMatrix:::.makeTempFile("folderOfLists")
dirMat  <- AnnotatedMatrix(dirPath)
dirMat$message("Loaded set of simple lists within a directory", bgcolor='cyan')

