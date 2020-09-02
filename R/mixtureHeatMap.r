
#' plot mixture of all two compartment profile combinations as panel; assumes eight compartments
#' @param Acup relative amount of a given cellular compartment protein that ends up in a given centrifugation fraction
#' @param totProt vector of amounts starting material in each fraction
#' @param eps small positive constant to add before taking a log transformation


mixtureHeatMap <- function(Acup=AcupMarkers, totProt, eps=0.01) {
  op <- par()

  par(mfrow=c(3,3))
  protAmtList <- c("Relative specific amount\n(RSA)",
                   "Normalized specific amount\n(NSA)",
                   "Relative amount\n(Acup)",
                   "Log2 RSA", "Log2 NSA", "Log2 Acup")
  protAmtMat <- matrix(protAmtList, nrow=2, byrow=T)
  library(grid)

  x <- c(0,1)
  y <- c(0,1)
  nComp <- nrow(Acup)  # number of compartments
  par(mfrow=c(nComp,nComp), mar=c(0.5, 0.5, 0.5, 0.5), xpd=NA)
 compartmentList <- rownames(Acup)
 # print column labels
 for (k in 1:nComp) {
   plot(y ~ x, type="n", axes=F, xlab="", ylab="")
   if (k != 1) text(0.5, 0.5, compartmentList[k], cex=2)
 }

 for (i in 1:(nComp-1)) {
   for (j in 1:nComp) {


   if (j <= i) {
     plot(y ~ x, type="n", axes=F, xlab="", ylab="")


     if ((j == 1) & (i != 8)) text(0.5, 0.5, compartmentList[i], cex=2)
   }
    if (j > i) {

  # i=1
  # j=4
  # create mixture
  mixProtiProtj <- proteinMix(AcupMarkers, Loc1=i, Loc2=j)
  mixProtiProtjRSA <- RSAfromAcup(Acup=mixProtiProtj$relAmount,
                                NstartMaterialFractions=6, totProt=tmtMS2totProt)

  # find RSA
  markerLocRrsa <- RSAfromS(SS=markerLocR, NstartMaterialFractions=6,
                          totProt=tmtMS2totProt)
  # apply CPA to relative specific amounts (recommended way)
  mixProtiProtjProp <- proLocAll(protProfileSummary=mixProtiProtjRSA,
                               markerLocR=markerLocRrsa, n.channels=9)

  # normalized specific amounts
  mixProtiProtjSpecAmt <- t(apply(mixProtiProtjRSA,1, function(x) x/sum(x)))

  # apply CPA to normalized specific amounts
  mixProtiProtjPropSpecAmt <- proLocAll(protProfileSummary=mixProtiProtjSpecAmt,
                                      markerLocR=markerLocR,
                                      n.channels=9)

  # apply CPA to relative amounts (Acup)
  mixProtiProtjPropAcup <-
    proLocAll(protProfileSummary=mixProtiProtj$relAmount,
            markerLocR=AcupMarkers, n.channels=9)

  # log transformed
  eps <- 0.00005


  mixProtiProtjPropLog2 <- proLocAll(protProfileSummary=log2(mixProtiProtjRSA + eps),
                                   markerLocR=log2(markerLocRrsa + eps), n.channels=9)


  mixProtiProtjPropSpecAmtLog2 <-
    proLocAll(protProfileSummary=log2(mixProtiProtjSpecAmt + eps),
            markerLocR=log2(markerLocR + eps), n.channels=9)

  mixProtiProtjPropAcupLog2 <-
    proLocAll(protProfileSummary=log2(mixProtiProtj$relAmount + eps),
            markerLocR=log2(AcupMarkers + eps) , n.channels=9)



# # # #

  ae11 <- mixtureAreaError(mixProtiProtjProp=mixProtiProtjProp, NstartMaterialFractions=6,
                 Loc1=i, Loc2=j, input.prop=mixProtiProtj$input.prop)
  ae12 <- mixtureAreaError(mixProtiProtjProp=mixProtiProtjPropSpecAmt,
                 NstartMaterialFractions=6, Loc1=i, Loc2=j,
                 input.prop=mixProtiProtj$input.prop)
  ae13 <- mixtureAreaError(mixProtiProtjProp=mixProtiProtjPropAcup, NstartMaterialFractions=6,
                 Loc1=i, Loc2=j, input.prop=mixProtiProtj$input.prop)



  ae21 <- mixtureAreaError(mixProtiProtjProp=mixProtiProtjPropLog2, NstartMaterialFractions=6,
                 Loc1=i, Loc2=j, input.prop=mixProtiProtj$input.prop)
  ae22 <- mixtureAreaError(mixProtiProtjProp=mixProtiProtjPropSpecAmtLog2, NstartMaterialFractions=6,
                 Loc1=i, Loc2=j, input.prop=mixProtiProtj$input.prop)
  ae23 <- mixtureAreaError(mixProtiProtjProp=mixProtiProtjPropAcupLog2,
                 NstartMaterialFractions=6, Loc1=i, Loc2=j,
                 input.prop=mixProtiProtj$input.prop)
  errorMat <- matrix(c(ae11, ae12, ae13, ae21, ae22, ae23), nrow=2, byrow=T)




  errorMat1m <- 2 - errorMat # so that higher errors are dark red
  library(plot.matrix)
  library("viridis")


  col <- rev(magma(20))
  plot(errorMat, col=col, breaks=seq(0, 2, 0.10), key=NULL, main="",
     axis.col=NULL, axis.row=NULL, xlab="", ylab="", digits=2, cex=0.8)
    }
   }
 }
 par(op)  # restore par

}


