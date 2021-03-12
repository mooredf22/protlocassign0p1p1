#' Compute area-based error for mixture of two compartment profiles
#' @param mixProtiProtjCPA data frame of CPA estimated proportions for each mixture
#' @param NstartMaterialFractions  Number of fractions that comprise the starting material
#' @param Loc1  row number of subcellular location 1 of mixture
#' @param Loc2  row number of subcellular location 2 of mixture

#' @param input.prop  true mixture proportions used to make the  mixtures (from proteinMix)
#' @param errorReturn  Return area of error region if true


mixtureAreaError <- function(mixProtiProtjCPA, NstartMaterialFractions=6, Loc1, Loc2,
                        input.prop) {
  # mixtures must be a list of equally spaced proportions
  # this program assumes exactly eight subcellular compartments
  # set up color and point lists
  Loc1 <- as.integer(Loc1)
  Loc2 <- as.integer(Loc2)
  loc.list <- names(mixProtiProtjCPA)

  fracList <- seq(0,1,0.1)

  areaErr <- 0
  for (kk in 1:ncol(mixProtiProtjCPA)) {
    #kk=1
    # use trapezoidal rule (trapz from package pracma)
    areaErr <- areaErr +
      abs(trapz(fracList, as.numeric(input.prop[,kk]) ) -
            trapz(fracList, as.numeric(mixProtiProtjCPA[,kk])))
  }

  return(areaErr)

}
