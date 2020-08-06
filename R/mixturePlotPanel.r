#' plot mixture of all two compartment profile combinations as panel; assumes eight compartments
#' @param Acup relative amount of a given cellular compartment protein that ends up in a given centrifugation fraction
#' @param totProt vector of amounts starting material in each fraction
#' @param fitType use (normalized) specific amounts or relative specific amounts
#' @param log2Transf use log2-transformed values
#' @param type normalized specific amounts or relative specific amounts
mixturePlotPanel <- function(Acup=AcupMarkers, totProt, fitType="relAmt", log2Transf=F) {
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
  # this program assumes exactly eight subcellular compartments
  # set up color and point lists
  loc.list <- names(mixProtiProtjProp)
  col.list <- c("red", "blue", "orange", "darkgreen", "orange", "lightblue", "purple", "green")
  pch.list <- c(1, 2, 3, 4, 17, 6, 15, 8)
  plotLables <- data.frame(loc.list, col.list, pch.list)
  par(mar=c(3,3,3,3))
  for (i in 1:7) {
    for (j in (i+1):8) {
      #mixturePlot(i=i, j=j, type=fitType, log2Transf=log2Transf)
      mixProtiProtj <- proteinMix(Acup=AcupMarkers, Loc1=i, Loc2=j)

      mixProtiProtjRSA <- RSAfromAcup(Acup=mixProtiProtj$mixAmount, NstartMaterialFractions=6,
                                      totProt=totProtUse)

      mixProtiProtjProp <- proLocAll(protProfileSummary=mixProtiProtjRSA,
                                     markerLocR=markerLocRrsa, n.channels=9)


      mixturePlot(mixProtiProtjProp=mixProtiProtjProp, NstartMaterialFractions=6, Loc1=i, Loc2=j,
                  mix.df=mixProtiProtj$mix.df)
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
}
