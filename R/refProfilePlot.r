

#' Plot profiles of reference genes
#'
#' This function plots profiles of reference genes and also the average profile for each compartment
#'
#' @param refLoc  the name of the reference subcellular compartment to plot
#' @param refLocProteins List of reference proteins
#' @param geneProfileSummary data frame of protein names and relative abundance levels.
#' @param markerLocR A matrix markerLocR giving the abundance level profiles of the subcellular locations
#'


refProfilePlot <- function(refLoc, refLocProteins=refLocProteinsJadot, geneProfileSummary=geneProfileSummaryTMTms2,
                           markerLocR=markerLocR)  {
  n.channels <- ncol(markerLocR)
  names(geneProfileSummary)[1] <- "geneName"
  meanReferenceGenes <- merge(x=refLocProteins, y=geneProfileSummary,
                            by.x="geneName", by.y="geneName", all.x=F, sort=F)
  markerLoc <- t(markerLocR)
  max.val <- max(markerLoc)
  #if (pdfout) pdf(file=paste("ClassificationMarkerProfilesTransformOutlierRej", dataUse, markersUse, ".pdf", sep=''), width=9, height=10)
  #if (!pdfout) windows(width=9, height=10)
  #if (dataUse == "TMT10revisit2") pdf(file="ClassificationMarkerProfiles2.pdf", width=9, height=10)
  #if (dataUse == "TMT10revisit3") pdf(file="ClassificationMarkerProfiles3.pdf", width=9, height=10)
  #par(mfrow=c(4,3))
  fractions.list <- rownames(markerLoc)
  compartments.list <- rownames(markerLocR)

  # stop if refLoc not found
  if (!(refLoc %in% compartments.list)) {
    cat("reference compartment not found\n")
    return(NULL)
  }

  location.list <- colnames(markerLoc)
  #location.list <- sort(unique(refLocProteins$referenceCompartment))
  #location.list2 <- unique(refLocProteins$referenceCompartment) # should be same as location.list
  #if (sum({location.list == location.list2}) != length(location.list)) {
  #    cat("error: locations in markerLoc don't match those in refLocProteins\n")
  #}


  #temp <- strsplit(fractions.list, "constrained")
  #unlist(temp)
  #fractions.list.short <- sub("constrained", "", fractions.list)


  n.loc <- ncol(markerLoc)
  #for (i in 1:n.loc) {
    # i=1
    #loc.i <- location.list[i]
    loc.i <- refLoc
    #assignLong.i <- assignmentLong.list[i]
    #channels.i <- meanCovarGenesStringent[meanCovarGenesStringent$AssignStringent == assign.i, 2+1:n.chan]
     channels.i <- meanReferenceGenes[meanReferenceGenes$referenceCompartment == loc.i,2+1:n.channels]
     refgenes.i <- as.character(meanReferenceGenes[meanReferenceGenes$referenceCompartment == loc.i,1])

     #channels.i <- nczfMeans[nczfMeans$AssignStringent == assign.i, c(3:8,10)]
    mean.i <- as.numeric(markerLoc[, {loc.i == location.list}] )
    xvals <- 1:length(mean.i)
    #if (!log2prop) mean.i <- 2^as.numeric(markerLoc[,i]) - eps
    #plot(mean.i ~ xvals, ylim=c(min.y,max.y.vec[i]), axes="F", type="l", ylab="")
    #plot(mean.i ~ xvals, ylim=c(-6,0.5), axes="F", type="l", ylab="")
    max.val <- max(channels.i)
    # ylim=c(0,max.val),
    plot(mean.i ~ xvals,ylim=c(0,max.val), axes="F", type="n", ylab="",
          xlab="")

    axis(1,at=xvals,labels=fractions.list, cex.axis=0.6)
    axis(2, las=1)
    for (j in 1:nrow(channels.i)) {
      means.j <- as.numeric(channels.i[j,])
      #if (!log2prop) means.j <- 2^as.numeric(channels.i[j,]) - eps
      lines(as.numeric(means.j) ~ xvals, col="red")
    }
    #lines(mean.i ~ xvals, lwd=1, lty=1, col="yellow")
    lines(mean.i ~ xvals, lwd=2, lty=2)
    title(paste(loc.i,"profiles"))
 # }




  }

