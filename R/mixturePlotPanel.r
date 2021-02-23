#' plot mixture of all two compartment profile combinations as panel; assumes eight compartments
#' Also assumes that "refLocationProfilesRSA" has been previously defined
#' @param Acup relative amount of a given cellular compartment protein that ends up in a given centrifugation fraction
#' @param totProt vector of amounts starting material in each fraction
#' @param fitType use (normalized) specific amounts or relative specific amounts
#' @param errorReturn return all area-based errors if true
#' @param log2Transf use log2-transformed values
#' @param type normalized specific amounts or relative specific amounts
mixturePlotPanel <- function(Acup=AcupMarkers, totProt, errorReturn=F, fitType="specAmt", log2Transf=F, eps=0.01) {
  #windows(width=8,height=10)
  layout(rbind(c(31,1,1,1,1),
               c(31,2,3,4,5),
               c(31,6,7,8,9),
               c(31,10,11,12,13),
               c(31,14,15,16,17),
               c(31,18,19,20,21),
               c(31,22,23,24,25),
               c(31,26,27,28,29),
               c(31,30,30,30,30)),

         heights=c(0.85,2,2,2,2,2,2, 2, .35),
         widths=c(0.4,2,2,2,2),respect=F)



  #layout.show(31)
  # this program assumes exactly eight subcellular compartments
  # set up color and point lists
  loc.list <- names(mixProtiProtjProp)
  col.list <- c("red", "blue", "orange", "darkgreen", "orange", "lightblue", "purple", "green")
  pch.list <- c(1, 2, 3, 4, 17, 6, 15, 8)

  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=F,cex=1)
  #aa <- 1.34
  #log2Transf <- F
  #log2Transf <- T
  #fitType <- "relAmt"
  #fitType <- "original"
  if (!log2Transf) transF <- "no transformation"
  if (log2Transf) transF <- "log2 transformed"
  text(x=2.5,y=0.4,paste("CPA Synthetic Protein Assignments", fitType, transF), cex=1.5)
  legend(x="bottom", legend=loc.list, col=col.list, pch=pch.list, ncol=8)

  # mixtures must be a list of equally spaced proportions

  plotLables <- data.frame(loc.list, col.list, pch.list)
  par(mar=c(3,3,3,3))
  kk <- 0
  mixErrorMat <- NULL
  for (i in 1:7) {
    for (j in (i+1):8) {
      #mixturePlot(i=i, j=j, type=fitType, log2Transf=log2Transf)
      mixProtiProtj <- proteinMix(Acup=AcupMarkers, Loc1=i, Loc2=j)

      mixProtiProtjRSA <- RSAfromAcup(Acup=mixProtiProtj$relAmount, NstartMaterialFractions=6,
                                      totProt=totProt)
      if (!log2Transf & {fitType == "rsa"}) {
       mixProtiProtjProp <- fitCPA(profile=mixProtiProtjRSA,
                                     refLocationProfiles=refLocationProfilesRSA, numDataCols=9)
      }

      if (log2Transf& {fitType == "rsa"}) {
      # Take a log2 transformation
        log2MixProtiProtjRSA <- log2(mixProtiProtjRSA + eps)
        log2MarkerLocRrsa <- log2(refLocationProfilesRSA + eps)
        mixProtiProtjProp <- fitCPA(profile=log2MixProtiProtjRSA,
                                     refLocationProfiles=log2MarkerLocRrsa, numDataCols=9)
      }
      if (!log2Transf & {fitType == "specAmt"}) {
        mixProtiProtjSpecific <- t(apply(mixProtiProtjRSA,1, function(x) x/sum(x)))
        markerLocRuse <- t(apply(markerLocRrsa,1, function(x) x/sum(x)))   # rows sum to one
        mixProtiProtjProp <- fitCPA(profile=mixProtiProtjSpecific,
                                       refLocationProfiles=markerLocRuse, numDataCols=9)
      }

      if (log2Transf& {fitType == "specAmt"}) {
        # Now convert to normalized specific amounts
        mixProtiProtjSpecific <- t(apply(mixProtiProtjRSA,1, function(x) x/sum(x)))
        markerLocRuse <- t(apply(refLocationProfilesRSA,1, function(x) x/sum(x)))  # rows sum to one
        # Take a log2 transformation
        log2MixProtiProtjSpecific <- log2(mixProtiProtjSpecific + eps)
        log2MarkerLocR <- log2(markerLocRuse + eps)


        mixProtiProtjProp <- fitCPA(profile=log2MixProtiProtjSpecific,
                                       refLocationProfiles=log2MarkerLocR, numDataCols=9)
      }

      if (!log2Transf & {fitType == "relAmt"}) {
        #mixProtiProtjSpecific <- t(apply(mixProtiProtjRSA,1, function(x) x/sum(x)))
        #markerLocRuse <- t(apply(markerLocRrsa,1, function(x) x/sum(x))) # rows sum to one
        mixProtiProtjProp <- fitCPA(profile=mixProtiProtj$relAmount,
                                       refLocationProfiles=AcupMarkers, numDataCols=9)
      }

      if (log2Transf& {fitType == "relAmt"}) {
        # Now convert to normalized specific amounts
        #mixProtiProtjSpecific <- t(apply(mixProtiProtjRSA,1, function(x) x/sum(x)))
        # Take a log2 transformation
        #log2MixProtiProtjSpecific <- log2(mixProtiProtj$relAmount + eps)
        #log2MarkerLocR <- log2(markerLocRuse + eps)


        mixProtiProtjProp <- fitCPA(profile=log2(mixProtiProtj$relAmount + eps),
                                       refLocationProfiles=log2(AcupMarkers + eps), numDataCols=9)
      }


      mixResult <- mixturePlot(mixProtiProtjProp=mixProtiProtjProp, NstartMaterialFractions=6, Loc1=i, Loc2=j,
                  input.prop=mixProtiProtj$input.prop, errorReturn=errorReturn)
      if (errorReturn) {
        mixErrorMat <- rbind(mixErrorMat, mixResult)
      }
    }
  }

  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=F,cex=1)
  #aa <- 1.34

  text(x=2.5,y=0.3,"True proportion", cex=2)

  x <- c(0,5)
  y <- c(0,0.5)
  par(mar=c(0,0,0,0))
  plot(y ~ x,type="n",axes=F,cex=1)
  #aa <- 1.34

  text(x=2.5,y=0.3,"Estimated proportion", cex=2, srt=90)
  if (errorReturn) return(mixErrorMat)
}
