#' plot mixture of two compartment profiles
#' @param Loc1 row number of the matrix Acup for compartment 1
#' @param Loc2 row number of the matrix Acup for compartment 2
#' @param mixProtiProtjProp data frame of CPA estimated proportions for each mixture
#' @param NstartMaterialFractions  Number of fractions that comprise the starting material
#' @param type normalized specific amounts or relative specific amounts

mixturePlot <- function(mixProtiProtjProp, NstartMaterialFractions=6, Loc1, Loc2, xaxisLab=T, yaxisLab=T) {
  # set up color and point lists
  loc.list <- names(mixProtiProtjProp)
  col.list <- c("red", "blue", "orange", "darkgreen", "orange", "lightblue", "purple", "green")
  pch.list <- c(1, 2, 3, 4, 17, 6, 15, 8)
  plotLables <- data.frame(loc.list, col.list, pch.list)

  fracList <- seq(0,1,0.1)  # for creating empty plot
  fracList2 <- fracList   #also need this for the empty plot

  # plot estimated proportion vs true proportion
  plot(fracList ~ fracList2, type="n", xlab="", ylab="", axes=F)  # this is the empty plot
  if (yaxisLab) axis(2, labels=T, at=c(0, 0.5, 1), las=1)
  if (!yaxisLab) axis(2, labels=F)
  if (xaxisLab) axis(1, labels=T, at=c(0, 0.5, 1))
  if (!xaxisLab) axis(1, labels=F)
  for (kk in 1:ncol(mixProtiProtjProp)) {
    points(mixProtiProtjProp[,kk] ~ fracList, col=col.list[kk], pch=pch.list[kk], cex=1.5)
  # calculate sum of squares of errors

  }
  area1 <- abs(0.5 + polyarea(x=c(fracList, 1), y=c(mixProtiProtjProp[,Loc1], 0)))
  area2 <- abs(-0.5 - polyarea(x=c(fracList, 0), y=c(mixProtiProtjProp[,Loc2],0)))
  area <- area1 + area2

  abline(a=0,b=1, col="gray")
  abline(a=1,b=-1, col="gray")
  titleText <- paste(loc.list[Loc1], "-", loc.list[Loc2])
  title(paste(titleText, "(", format(round(area,3), nsmall=3), ")"))   # guarantee 3 digits after decimal
  #plotLables
  #loc.list  col.list pch.list
  #1     Cyto       red        1
  #2       ER      blue        2
  #3    Golgi    orange        3
  #4     Lyso darkgreen        4
  #5     Mito    orange       17
  #6      Nuc lightblue        6
  #7    Perox    purple       15
  #8       PM     green        8
}
