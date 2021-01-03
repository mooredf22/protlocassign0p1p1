#' Plot profiles of reference proteins
#'
#' This function plots the average profiles of any protein in the dataset,
#'   the peptide profiles, and also the reference profile for each compartment
#'
#' @param protName Name of the protein to plot
#' @param refLocProteins List of reference proteins
#' @param protProfileSummary data frame of protein names and their relative abundance levels..
#' @param Nspectra indicator for if there are columns in protProfileSummary for Nspectra
#'         (number of spectra) and Npep (number of peptides)
#' @param finalList spectrum-level abundance levels by protein and peptide; Ehis is NULL
#'        if not available
#' @param n.fractions  number of fractions per protein
#' @param n.compartments number of compartments (8 in Jadot data)
#' @param markerLocR A matrix markerLocR giving the abundance level profiles of the subcellular locations
#'        n.compartments = 8 columns are subcellular locations, and n.fractions rows are the fraction names
#' @param assignPropsMat A matrix of assignment proportions, from the constrained proportional assignment algorithm,
#'        and optionally upper and lower 95 percent confidence limits
#' @param propCI True if lower and upper confidence intervals are included in assignPros
#' @param predictedProp.mat A matrix of CPA predicted proportions, from proLocAll
#' @param confint indicator for if there are standard errors (from bootstrapping) for the assigned proportions
#'        (in protProfileSummary)
#' @param predictedPropL.mat Matrix of lower confidence limits, one row for each protein, and 9 columns
#'        of lower confidence interval limits
#' @param predictedPropU.mat Matrix of upper confidence limits, one row for each protein, and 9 columns
#'        of upper confidence interval limits
#'

protPepPlotfun <- function(protName, protProfileSummary, Nspectra=T, finalList=NULL,
                        protPepData=NULL, refCol=5, n.fractions=9, n.compartments=8,
                        markerLocR=markerLocRuse, assignPropsMat=assignPropsUse, propCI=F,
                        transType="") {
  # protPlot is the number of the protein to plot
  # protProfileSummaryUse is a matrix with components:
  #    rownames: name of protein
  #    N, M, L1, ... : relative protein levels for fractions 1 through n.fractions; must sum to 1
  #
  # If standard errors are not available, se is false, and seMat is NULL
  #  If standard errors are available, se is true, and seMat is the list of proteinss and n.fraction standard errors
  # If Nspectra and Npep (number of peptides) are included, Nspectra=T
  #  markerLoc:  matrix with n.fractions rows and 8 columns.
  #      Column names are the subcellular fractions, Cytosol, ER, Golgi, etc.
  #      Row names are the names of the fractions: N, M, L1, L2, etc.
  #      Column and row names are required
  # protPlot <- 1
  # protPlot <- ind.prot
  # protPlot <- 6540
  # protPlot <- 93
  #protName.i <- as.character(meanCovarProtsUseAll.t$protNames[protPlot])
  # finalList contains all peptides AND spectra
  # protPepData contains protein profiles and mean peptide profiles;
  #   If this is provided, you cannot have either finalList or protProfileSummary
  # Can't have both!
  #if (!is.null(finalList) & !is.null(protPepData)) {
  #  cat("Cannot have both finalList and protPepData\n")
  #  stop()
  #}
  #if (!is.null(protProfileSummary) & !is.null(protPepData)) {
  #  cat("Cannot have both protProfileSummary and protPepData\n")
  #  stop()
  #}
  #if (!is.null(protPepData)) {
    # extract protProfileSumary and pepProfileSummary from protPepData
   # protPepData.i <- protPepData[protPepData$prot == protName,]
  #}
  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  protsOK <- {rownames(protProfileSummary) == rownames(assignPropsMat)}

  temp <- protIndex(protName, protProfileSummary, exactMatch=T)
  protName.i <- protName
  # must be exact to avoid duplicate finds

  # this can be a vector, matrix, or vector, so handling is complicated
  #if (is.matrix(tempx)) temp <- tempx   # leave it alone
  # if temp is not a matrix, can then test for being NA with no error returned
  if (!is.matrix(temp)) {
    if (is.na(temp)[1]) {
      cat(paste(protName, " not found \n"))  # first element is Na
      return(temp)
    }
  }
  if (nrow(temp) > 1) {
    cat(paste("more than one protein matches pattern \n"))
    return(temp)
  }
  protPlot <- temp[1,1]   # works even if temp is a vector
  # protPlot is the index number of the protein

  ##if(!protsOK) cat("Error: protein names don't match\n")
  ##stopifnot(protsOK)
  #assignProbsOut <- protProfileSummary[,1:(1+8)]   # just the protein name and the assigned proportions to the 8 compartments
  meanProteinLevels <- protProfileSummary[,1:n.fractions]  # just the profiles

  protNames <- rownames(protProfileSummary)
  protNameUnique <- make.unique(protNames)
  rownames(meanProteinLevels) <- protNameUnique
  subCellNames <- rownames(markerLocR)

  fractions.list <- colnames(markerLocR)  # column names[protPlot])

  # # # # # # # # # # # # # #
  #  Do the following if "finalList" (the full list of peptides and spectra)
  #   are available
  # # # # # # # # # # # # # #

  #if (!is.null(finalList)) {


  # # # # # # # # # # # # # #
  #  Do the following if "peptideList" (the list of mean peptide profiles)
  #   is available
  # # # # # # # # # # # # # #
  if (!is.null(protPepData)){

    protPepData.i <- protPepData[protPepData$prot == protName.i,]
    means.peptides.i <- protPepData.i[,refCol + 1:9]
    n.spectra.i <- protPepData.i$Nspectra.1
    outlierFlagVec.i <- {protPepData.i$outlier.num.peptides > 0}
    n.unique.peptide.i <- nrow(protPepData.i)
  }


  yy <- as.numeric(meanProteinLevels[protPlot,])

  xvals <- 1:n.fractions
  #max.y <- max(channelsAll, na.rm=T)
  #max.y <- max(meanProteinLevels, na.rm=T)
  min.y <-0
  #if (dataType == "raw") min.y <- 0
  loc.list <- subCellNames
  #windows(width=10, height=8)
  #par(mfrow=c(3,3))
  layout(rbind(c(14,1,1,1,15),
               c(14,2,3,4,15),
               c(14,5,5,5,15),
               c(14,6,7,8,15),
               c(14,9,9,9,15),
               c(14,10,11,12,15),
               c(14,13,13,13,15)),
         #c(16,12,13,14,15,17),
         #c(16,18,18,18,18,17)),
         heights=c(0.75,2,0.25,2,0.25,2,0.25),
         widths=c(0.4,2,2,2,0.4),respect=F)



  #layout.show(15)
  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=F,cex=1)
  #aa <- 1.34

  #if (length(indAssignProp.prot) > 0) {
    text(x=2.5,y=0.3,paste(protName.i), cex=2)

    NpeptidesPlot <- protProfileSummary$Npep[protPlot]
    NspectraPlot <- protProfileSummary$Nspectra[protPlot]
    NpeptidesPlotText <- " peptides and "

    #browser()

    if (NpeptidesPlot == 1) NpeptidesPlotText <- " peptide and "
    NspectraPlotText <- " spectra"
    if (NspectraPlot == 1) NspectraPlotText <- " spectrum"
    text(x=2.5,y=0.1, paste(NpeptidesPlot, NpeptidesPlotText,
                            NspectraPlot , NspectraPlotText), cex=2)
  #}
 # if (length(indAssignProp.prot) ==0) {
 #   text(x=2.5,y=0.3,protName.i, cex=2)
 # }

  # max.y <- max(c(max(means.peptides.i), max(markerLocR[i,])))
  min.y <- 0
  par(mar=c(2,4,2,1.5))
  # The plots are in alphabetical order
  # re-arrange the plots so that they are in this order:
  #   Mito (7)     Lyso  (4)   Perox (5)
  #   ER (2)       Golgi (1)   PM (8)
  #   Cyto (3 )    Nuc (6)
  #
  #loc.ord <- c(7, 4, 5, 2, 1, 8, 3, 6)
  loc.ord <- c(5, 4, 7, 2, 3, 8, 1, 6)
  for (i in loc.ord) {    # do all the subcellular locations
    # i=1
    if (T) {
      #if ({loc.ord[i] == 4} | {loc.ord[i] == 7}) {
      if ({i == 2} | {i == 1}) {
        x <- c(0,5)
        y <- c(0,0.5)
        par(mar=c(0,1,0,0))
        plot(y ~ x,type="n",axes=F,cex=1, ylab="")
        par(mar=c(2,3.5,2,1.5))
      }
    }

    assign.i <- names(meanProteinLevels)[i]  # channel i name

    #assignLong.i <- names(markerLoc)[i]
    assignLong.i <- subCellNames[i]
    #channels.i <- meanProteinLevels[{names(as.data.frame(markerLocR)) == as.character(assign.i)},]
    ##?? means.peptides.i <- meanProteinLevels[{names(as.data.frame(markerLocR)) == as.character(assign.i)},]
    #mean.i <- markerLoc[,i]
    mean.i <- as.numeric(markerLocR[i,])
    max.y <- max(c(max(means.peptides.i), max(markerLocR[i,])))
    n.uniq.peptide.i <- nrow(means.peptides.i)
    #if (!is.null(finalList)) max.y <- max(c(max(means.peptides.i), max(markerLocR[i,])))
    #if (is.null(finalList)) max.y <- max(c(mean.i,yy))
    par(mar=c(2,4.1,2,1.5))
    plot(mean.i ~ xvals,  axes="F", type="l",
         ylim=c(min.y, max.y), ylab=transType)
    axis(1,at=xvals,labels=fractions.list)
    axis(2)

  if (T) {
    for (j in 1:n.uniq.peptide.i) {
      lwdplot <- 1
      colplot <- "cyan"
      if (n.spectra.i[j] > 1) {
        lwdplot <- 2
        colplot <- "deepskyblue"
      }
      if (n.spectra.i[j] > 2) {
        lwdplot <- 3
        colplot <- "dodgerblue3"
      }
      if (n.spectra.i[j] > 5)   {
        lwdplot <- 4
        colplot <- "blue"
      }

      lines(as.numeric(means.peptides.i[j,]) ~ xvals, cex=0.5, lwd=lwdplot, col=colplot)
      if (outlierFlagVec.i[j] == 1) lines(as.numeric(means.peptides.i[j,]) ~ xvals,
                                          cex=0.5, lwd=2, col="orange")

    }
  }




    lines(yy ~ xvals, col="red", lwd=2)
    lines(mean.i ~ xvals, lwd=2, col="yellow")
    lines(mean.i ~ xvals, lwd=2, col="black", lty=2)
    if (propCI) {
      predPrL.i <- predictedPropL.mat[indAssignProb.prot,i]
      predPrU.i <- predictedPropU.mat[indAssignProb.prot,i]

      #  if all standard errors are zero, don't make confidence limits
      allZero <- {sum(seMat.i, na.rm=T) == 0}
      if (allZero) {
        predPrL.i <- NA
        predPrU.i <- NA
      }

      title(paste(assignLong.i, "\n p = ",
                  round(channelsAll[indAssignProb.prot,i], digits=2 ),
                  " (",
                  round(predPrL.i, digits=2 ),
                  ", ",
                  round(predPrU.i, digits=2 ),
                  ")"
      ))
    }
    if (!propCI) {
      title(paste(assignLong.i, "\n p = ", round(assignPropsMat[protPlot,i], digits=2 )))
    }
  }
  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x, type="n", axes=F)
  if (T) {
    legend(x=1, y=0.4, legend=c("Reference profile", "Average profile", "1 spectrum", "2 spectra", "3-5 spectra", "6+ spectra", "Outliers"),
         col=c("yellow", "red", "cyan", "deepskyblue", "dodgerblue3", "blue", "orange"),
         lwd=c(2,2,1,2,3,4,2), lty=c(1,1,1,1,1,1,1))
    legend(x=1, y=0.4, legend=c("Reference profile", "Average profile", "1 spectrum", "2 spectra", "3-5 spectra", "6+ spectra", "Outliers"),
         col=c("black", "red", "cyan", "deepskyblue", "dodgerblue3", "blue", "orange"),
         lwd=c(2,2,1,2,3,4,2), lty=c(2,1,1,1,1,1,7))
  }
  if (F) {
    legend(x=1, y=0.4, legend=c("Reference profile", "Average profile"),
           col=c("yellow", "red"), lwd=c(2,2), lty=c(1,1))
    legend(x=1, y=0.4, legend=c("Reference profile", "Average profile"),
           col=c("black", "red"), lwd=c(2,2), lty=c(2,1))
  }
}
