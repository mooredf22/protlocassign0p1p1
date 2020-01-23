#' Plot profiles of reference genes
#'
#' This function plots the average profiles of any gene in the dataset,
#'   the peptide profiles, and also the reference profile for each compartment
#'
#' @param refLocProteins List of reference proteins
#' @param geneProfileSummary data frame of protein names and relative abundance levels.
#' @param finalList spectrum-level abundance levels by protein and peptide
#'
#' @param n.fractions  number of fractions per protein
#' @param Nspectra indicator for if there are columns in geneProfileSummary for Nspectra
#'         (number of spectra) and Nseq (number of peptides)
#' @param matLocR A matrix matLocR giving the abundance level profiles of the subcellular locations
#'        8 columns are subcellular locations, and n.fractions rows are the fraction names
#' @param confint indicator for if there are standard errors (from bootstrapping) for the assigned proportions
#'        (in geneProfileSummary)
#' @param predictedPropL.mat Matrix of lower confidence limits, one row for each gene, and 9 columns
#'        of lower confidence interval limits
#' @param predictedPropU.mat Matrix of upper confidence limits, one row for each gene, and 9 columns
#'        of upper confidence interval limits
#'

protPlotfun <- function(protPlot, geneProfileSummary=geneProfileSummaryUse, finalList=finalListUse, n.fractions=9,
                        Nspectra=T, matLocR=matLocRuse, confInt=F, predictedPropL.mat=NULL, predictedPropU.mat) {
  # protPlot is the number of the protein (gene) to plot
  # geneProfileSummaryUse is a matrix with components:
  #    geneName: name of gene (protein)
  #    N, M, ... : assignment proportions for fractions 1 through n.fractions
  #
  # If standard errors are not available, se is false, and seMat is NULL
  #  If standard errors are available, se is true, and seMat is the list of genes and n.fraction standard errors
  # If Nspectra and Nseq (number of peptides) are included, Nspectra=T
  #  matLoc:  matrix with n.fractions rows and 8 columns.
  #      Column names are the subcellular fractions, Cytosol, ER, Golgi, etc.
  #      Row names are the names of the fractions: N, M, L1, L2, etc.
  #      Column and row names are required
  # protPlot <- 1
  # protPlot <- ind.prot
  # protPlot <- 6540
  # protPlot <- 93
  #protName.i <- as.character(meanCovarGenesUseAll.t$protNames[protPlot])
  assignProbsOut <- geneProfileSummary[,1:(1+n.fractions)]   # just the gene name and the n.fractions assigned proportions
  meanCovarGenes <- assignProbsOut[,-1]  # drop gene name column
  channelsAll <-  assignProbsOut[,-1]  # another name for this
  row.names(meanCovarGenes) <- assignProbsOut$geneName
  subCellNames <- names(matLoc)
  subCellNames[8] <- "Plasma membrane"
  subCellNames[4] <- "Lysosomes"
  subCellNames[2] <- "Endoplasmic reticulum"

  matLoc <- t(matLocR)
  fractions.list <- row.names(matLoc)

  protName.i <- as.character(assignProbsOut$geneName[protPlot])
  if (se) seMat.i <- seMat[protPlot,]  # just the standard errors
  #finalList.t <- finalList
  finalList.i <- finalList[toupper(finalList$gene) == protName.i, ]
  fractions.i <- finalList.i[,2 + c(1:n.fractions)]
  #  add these back in:
  #exc.ind.i <- {finalList.i$outCountBoxPlot > 0}
  #fractions.use.i <- fractions.i[!exc.ind.i,]
  #finalList.use.i <- finalList.i[!exc.ind.i,]
  finalList.use.i <-  finalList.i
  fractions.use.i <- fractions.i
  peptide.i <- as.character(finalList.use.i$peptide)
  n.uniq.peptide.i <- length(unique(peptide.i))
  uniq.peptides.list <- unique(peptide.i)
  means.peptides.i <- matrix(NA, nrow=n.uniq.peptide.i, ncol=n.fractions)
  n.spectra.i <- rep(NA, n.uniq.peptide.i)
  #browser()
  for (jj in 1:n.uniq.peptide.i) {
    fractions.use.i.jj <- fractions.use.i[uniq.peptides.list[jj] == peptide.i,]
    means.peptides.i[jj,] <- apply(fractions.use.i.jj,2,mean)
    n.spectra.i[jj] <- nrow(fractions.use.i.jj)
  }

  n.assign <- nrow(assignProbsOut)
  indAssignProb.prot <- (1:n.assign)[assignProbsOut$geneName == protName.i]


  n.fractions.i <- nrow(fractions.i)


  #yy= as.numeric(channelsAll[protPlot,])
  yy <- as.numeric(meanCovarGenes[protPlot,])

  xvals <- 1:n.fractions
  #max.y <- max(channelsAll, na.rm=T)
  max.y <- max(meanCovarGenes, na.rm=T)
  min.y <-0
  #if (dataType == "raw") min.y <- 0
  loc.list <- names(matLoc)
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
  if (length(indAssignProb.prot) > 0) {
    text(x=2.5,y=0.3,paste(protName.i), cex=2)

    NpeptidesPlot <- geneProfileSummary$Nseq[indAssignProb.prot]
    NspectraPlot <- geneProfileSummary$Nspectra[indAssignProb.prot]
    NpeptidesPlotText <- " peptides and "
    if (NpeptidesPlot == 1) NpeptidesPlotText <- " peptide and "
    NspectraPlotText <- " spectra"
    if (NspectraPlot == 1) NspectraPlotText <- " spectrum"
    text(x=2.5,y=0.1, paste(NpeptidesPlot, NpeptidesPlotText,
                            NspectraPlot , NspectraPlotText), cex=2)
  }
  if (length(indAssignProb.prot) ==0) {
    text(x=2.5,y=0.3,protName.i, cex=2)
  }

  par(mar=c(2,1.5,2,1.5))
  for (i in 1:8) {    # do all the subcellular locations
    # i=1
    if (T) {
      if ({i == 4} | {i == 7}) {
        x <- c(0,5)
        y <- c(0,0.5)
        par(mar=c(0,0,0,0))
        plot(y ~ x,type="n",axes=F,cex=1)
        par(mar=c(2,1.5,2,1.5))
      }
    }

    assign.i <- names(channelsAll)[i]  # channel i name

    #assignLong.i <- names(matLoc)[i]
    assignLong.i <- subCellNames[i]
    channels.i <- channelsAll[{names(as.data.frame(t(matLoc))) == as.character(assign.i)},]
    mean.i <- matLoc[,i]
    plot(mean.i ~ xvals,  axes="F", type="l",
         ylim=c(min.y, max.y))
    axis(1,at=xvals,labels=fractions.list)
    axis(2)

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

    }




    lines(yy ~ xvals, col="red", lwd=2)
    lines(mean.i ~ xvals, lwd=2, col="yellow")
    lines(mean.i ~ xvals, lwd=2, col="black", lty=2)
    if (se) {
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
    if (!se) {
      title(paste(assignLong.i, "\n p = ", round(channelsAll[indAssignProb.prot,i], digits=2 )))
    }
  }
  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x, type="n", axes=F)
  legend(x=1, y=0.4, legend=c("Reference profile", "Average profile", "1 spectrum", "2 spectra", "3-5 spectra", "6+ spectra"),
         col=c("yellow", "red", "cyan", "deepskyblue", "dodgerblue3", "blue"), lwd=c(2,2,1,2,3,4), lty=c(1,1,1,1,1,1))
  legend(x=1, y=0.4, legend=c("Reference profile", "Average profile", "1 spectrum", "2 spectra", "3-5 spectra", "6+ spectra"),
         col=c("black", "red", "cyan", "deepskyblue", "dodgerblue3", "blue"), lwd=c(2,2,1,2,3,4), lty=c(2,1,1,1,1,1))
}
