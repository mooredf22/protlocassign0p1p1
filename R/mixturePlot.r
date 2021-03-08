#' plot mixture of two compartment profiles
#' @param mixProtiProtjCPA data frame of CPA estimated proportions for each mixture
#' @param NstartMaterialFractions  Number of fractions that comprise the starting material
#' @param Loc1  row number of subcellular location 1 of mixture
#' @param Loc2  row number of subcellular location 2 of mixture

#' @param input.prop  true mixture proportions used to make the  mixtures (from proteinMix)
#' @param errorReturn  Return area of error region if true

mixturePlot <- function(mixProtiProtjCPA, NstartMaterialFractions=6, Loc1, Loc2,
                        input.prop, errorReturn=F, subTitle=NULL, xaxisLab=T, yaxisLab=T) {
  # mixtures must be a list of equally spaced proportions
  # this program assumes exactly eight subcellular compartments
  # set up color and point lists
  Loc1 <- as.integer(Loc1)
  Loc2 <- as.integer(Loc2)
  loc.list <- names(mixProtiProtjCPA)
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
  for (kk in 1:ncol(mixProtiProtjCPA)) {
    points(mixProtiProtjCPA[,kk] ~ fracList, col=col.list[kk], pch=pch.list[kk], cex=0.85)
  # calculate sum of squares of errors

  }

  areaErr <- 0
  for (kk in 1:ncol(mixProtiProtjCPA)) {
    #kk=1
    # use trapezoidal rule (trapz from package pracma)
    areaErr <- areaErr +
                abs(trapz(fracList, as.numeric(input.prop[,kk]) ) -
                      trapz(fracList, as.numeric(mixProtiProtjCPA[,kk])))



  }

  segments(0,0,1,1, col=col.list[Loc1])
  segments(0,1,1,0, col=col.list[Loc2])
  segments(0,0,1,0, col="gray")
  titleText <- paste(loc.list[Loc1], "-", loc.list[Loc2])
  title(paste(titleText, "(", format(round(areaErr,3), nsmall=3), ")\n", subTitle))   # guarantee 3 digits after decimal
  #plotLables
  #loc.list  col.list pch.list
  #1     Cyto       red        1  open circle
  #2       ER      blue        2  triangle
  #3    Golgi    orange        3  plus
  #4     Lyso darkgreen        4  X
  #5     Mito    orange       17  solid triangle
  #6      Nuc lightblue        6  upside down triangle
  #7    Perox    purple       15  solid square
  #8       PM     green        8  asterisk
  if (errorReturn) {
    areaErrOut <- data.frame(loc.list[Loc1], loc.list[Loc2], areaErr)
    names(areaErrOut) <- c("Loc1", "Loc2", "ErrorArea")
    return(areaErrOut)
  }
}
