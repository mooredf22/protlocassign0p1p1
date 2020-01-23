

#' Plot profiles of reference genes
#'
#' This function plots profiles of reference genes and also the average profile for each compartment
#'
#' @param refLocProteins List of reference proteins
#' @param geneProfileSummary data frame of protein names and relative abundance levels.
#' @param matLocR A matrix matLocR giving the abundance level profiles of the subcellular locations
#'


referenceProfilePlot <- function(refLocProteins=refLocProteinsJadot, geneProfileSummary=geneProfileSummaryTMTms2, matLocR=matLocR,
                                 n.channels=7, dataUse="", markersUse="")  {
  names(geneProfileSummary)[1] <- "geneName"
  meanReferenceGenes <- merge(x=refLocProteins, y=geneProfileSummary,
                            by.x="geneName", by.y="geneName", all.x=F, sort=F)
  matLoc <- t(matLocR)
  max.val <- max(matLoc)
  #if (pdfout) pdf(file=paste("ClassificationMarkerProfilesTransformOutlierRej", dataUse, markersUse, ".pdf", sep=''), width=9, height=10)
  #if (!pdfout) windows(width=9, height=10)
  #if (dataUse == "TMT10revisit2") pdf(file="ClassificationMarkerProfiles2.pdf", width=9, height=10)
  #if (dataUse == "TMT10revisit3") pdf(file="ClassificationMarkerProfiles3.pdf", width=9, height=10)
  #par(mfrow=c(4,3))
  fractions.list <- rownames(matLoc)

  #location.list <- colnames(matLoc)
  location.list <- sort(unique(refLocProteins$referenceCompartment))


  #temp <- strsplit(fractions.list, "constrained")
  #unlist(temp)
  fractions.list.short <- sub("constrained", "", fractions.list)


  n.loc <- ncol(matLoc)
  for (i in 1:n.loc) {
    # i=1
    loc.i <- location.list[i]
    #assignLong.i <- assignmentLong.list[i]
    #channels.i <- meanCovarGenesStringent[meanCovarGenesStringent$AssignStringent == assign.i, 2+1:n.chan]
     channels.i <- meanReferenceGenes[meanReferenceGenes$referenceCompartment == loc.i,2+1:n.channels]

    #channels.i <- nczfMeans[nczfMeans$AssignStringent == assign.i, c(3:8,10)]
    mean.i <- as.numeric(matLoc[,i] )
    xvals <- 1:length(mean.i)
    #if (!log2prop) mean.i <- 2^as.numeric(matLoc[,i]) - eps
    #plot(mean.i ~ xvals, ylim=c(min.y,max.y.vec[i]), axes="F", type="l", ylab="")
    #plot(mean.i ~ xvals, ylim=c(-6,0.5), axes="F", type="l", ylab="")
    # max.val <- 0.7
    plot(mean.i ~ xvals, ylim=c(0,max.val), axes="F", type="n", ylab="",
          xlab="")
    axis(1,at=xvals,labels=fractions.list, cex.axis=0.6)
    axis(2, las=1)
    for (j in 1:nrow(channels.i)) {
      means.j <- as.numeric(channels.i[j,])
      #if (!log2prop) means.j <- 2^as.numeric(channels.i[j,]) - eps
      lines(as.numeric(means.j) ~ xvals, col="red")
    }
    lines(mean.i ~ xvals, lwd=3, lty=1, col="yellow")
    lines(mean.i ~ xvals, lwd=3, lty=2)
    title(paste(loc.i,"profiles"))
  }

#  if (pdfout) dev.off()
  }

