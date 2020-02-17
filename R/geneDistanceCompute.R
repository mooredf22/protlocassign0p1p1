# compute distances, and find nearest genes

# first, create a distance matrix, n.genes by n.genes
#distUse <- dist(geneProfileSummaryTMTms2[,2:10], method="euclidean")
#genesUse <- geneProfileSummaryTMTms2[,1]
#geneName <- "AADAC"

#' Compute distances between a particular gene or profile and all other genes, and list the nearest ones
#'
#' @param geneName  Name of gene to which distances are to be computed
#' @param n.nearest Number of nearest genes to list
#' @param distGenes distance matrix created by, for example, dist(geneProfileSummaryTMTms2[,2:10])
#' @param geneNames A list of all genes in a dataset such as geneProfileSummaryTMTms2[,1]
#' @return List of the genes in geneName closest to geneName or to the profile
#' @examples
#' distUse <- dist(geneProfileSummaryTMTms2[,2:10], method="euclidean")
#' genesUse <- geneProfileSummaryTMTms2[,1]
#' nearestGenes("AADAC", n.nearest=10,  distGenes=distUse, geneNames=genesUse)
nearestGenes <- function(geneName, n.nearest=5, distGenes=distUse, geneNames=genesUse) {
  distGenesMat <- as.matrix(distGenes)

    ref <- protIndex(geneName)
    if (is.character(ref)) {
      return(ref)
    }
    if (nrow(ref) > 1) {
      cat("More than one gene matches geneName\n")
      return(ref)
    }

    ind.ref <- ref[1,1]
    vect.dist <- distGenesMat[ind.ref,]  # vector of distances to the reference gene



  nearest.list <- sort(vect.dist)
  nearest.indices <- as.numeric(names(nearest.list))[1:n.nearest]
  result <- data.frame(geneNames[nearest.indices], vect.dist[nearest.indices])
  names(result) <- c("geneName", "euclidean distance")
  result

}

#nearestGenes("AADAC", n.nearest=10,  distGenes=distUse, geneNames=genesUse)

#nearestGenes("AAD", n.nearest=10,  distGenes=distUse, geneNames=genesUse)

