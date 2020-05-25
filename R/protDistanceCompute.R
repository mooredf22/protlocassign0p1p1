# compute distances, and find nearest genes

# first, create a distance matrix, n.genes by n.genes
#distUse <- dist(protProfileSummaryTMTms2[,2:10], method="euclidean")
#protsUse <- protProfileSummaryTMTms2[,1]
#protName <- "AADAC"

#' Compute distances between a particular gene or profile and all other genes, and list the nearest ones
#'
#' @param protName  Name of protein to which distances are to be computed
#' @param n.nearest Number of nearest proteins to list
#' @param distProts distance matrix created by, for example, dist(protProfileSummaryTMTms2[,2:10])
#' @param protNames A list of all proteins in a dataset such as protProfileSummaryTMTms2[,1]
#' @return List of the proteins in protName closest to protName or to the profile
#' @examples
#' distUse <- dist(protProfileSummaryTMTms2[,2:10], method="euclidean")
#' protsUse <- protProfileSummaryTMTms2[,1]
#' nearestProts("AADAC", n.nearest=10,  distProts=distUse, protNames=protsUse)
nearestProts <- function(protName, n.nearest=5, distProts=distUse, protNames=protsUse) {
  distProtsMat <- as.matrix(distProts)

    ref <- protIndex(protName)
    if (is.character(ref)) {
      return(ref)
    }
    if (nrow(ref) > 1) {
      cat("More than one protein matches protName\n")
      return(ref)
    }

    ind.ref <- ref[1,1]
    vect.dist <- distProtsMat[ind.ref,]  # vector of distances to the reference protein



  nearest.list <- sort(vect.dist)
  nearest.indices <- as.numeric(names(nearest.list))[1:n.nearest]
  result <- data.frame(protNames[nearest.indices], vect.dist[nearest.indices])
  names(result) <- c("protName", "euclidean distance")
  result

}

#nearestProts("AADAC", n.nearest=10,  distProts=distUse, protNames=protsUse)

#nearestProts("AAD", n.nearest=10,  distProts=distUse, protNames=genesUse)

