#' Set Fisher Example Object
#'
#' Create a simple toy SetFisher object
#'
#' @details
#'
#' The object created can be used to explore SetFisher
#' functionality. It is also used by the test framework for the
#' package.
#'
#' The object will have two SetFisherAnalysis objects associated with
#' it. Both are silly fictional examples. The Query matrix for both
#' analyses is a set of cyber attacks by various super villains.
#'
#' In one analysis the Ontology is directly associated with the
#' villain code names (so no Mapping matrix is required), and
#' associates villains with various corporate endorsements. This
#' analysis could be used to see if any particular attack has an
#' unusual representation of villains paid by a corporation.
#'
#' The other analysis is using an ontology that utilizes Villain
#' Identifier Numbers (VIN- accessions). This is a more robust method
#' of tracking villains, since they often use the same code names. It
#' requires a Mapping matrix however, which maps the code names to VIN
#' accessions. This then allows use of the Criminal Organization
#' ontology, which can be used to see if any organization was
#' unusually represented at a particular crime.
#'
#' The matrices are intentionally small in order to aid manual
#' inspection. Because of this the enrichment p-values are going to be
#' insignificant by normal criteria.
#'
#' @param usemap Default \code{TRUE} which will return a SetFisher
#'     object that hs (requires) the optional Mapping Matrix in order
#'     to connect the Query Matrix to the Ontology Matrix
#'
#' @examples
#'
#' esf <- setFisherExampleObject()
#'
#' @importFrom AnnotatedMatrix annotatedMatrixExampleFile
#' 
#' @export

setFisherExampleObject <- function(usemap=TRUE) {
    sf <- SetFisher()
    crimes <- AnnotatedMatrix::annotatedMatrixExampleFile(
        'CrimeScenes.mtx', pkg='SetFisher')
    sf$defaultQuery( crimes )
    if (usemap) {
        name2vin <- AnnotatedMatrix::annotatedMatrixExampleFile(
            'VillainCodeNames.mtx', pkg='SetFisher')
        orgs <- AnnotatedMatrix::annotatedMatrixExampleFile(
            'CriminalOrganizations.mtx', pkg='SetFisher')
        sf$analysis(orgs, idmap=name2vin)
    } else {
        endorse <- AnnotatedMatrix::annotatedMatrixExampleFile(
            'VillainEndorsements.mtx', pkg='SetFisher')
        sf$analysis(endorse)
    }
    sf
}
