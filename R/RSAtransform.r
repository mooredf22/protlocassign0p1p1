

#' Set up RSA (relative specific amount) profiles for constrained proportional assignment
#'
#' @param NSA A matrix giving the specific amount or normalized specific amount
#'           profiles, either marker profiles from 'cpaSetup' (which requires
#'           normalized specific amounts as input), or a list of many protein profiles
#' @param NstartMaterialFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param totProt Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's
#'
#' @return amtProtFrac: amount of protein in a given fraction
#' @return relAmtProtFrac amount of given protein in fraction / amount of given protein in starting material

# # # # # # # # #
# name changes July 8, 2020
# # # # # # # # #

# function name change
#abundanceTransform <- function(markerProfiles, NstartMaterialFractions=6,

# markerProfiles -> NSA   # no need to limit this to the marker profiles
# amtProtFrac -> AA  # these are the amounts of protein (or other species)
# relAmtProtFrac -> Acup  # relative amounts
                # normalize specific amounts in each fraction
                # to the amount in starting material
    # these are of interest for simulations
#

#relAmtTransform <- function(NSA, NstartMaterialFractions=6,
AcupFromNSA <- function(NSA, NstartMaterialFractions=6,
                               totProt=NULL) {
  if (ncol(NSA) != length(totProt)) {
    cat("Error from relAmtTransoform: no. of rows of NSA must match length of totProt\n")
  }

  startMaterialFractions <- NSA[,1:NstartMaterialFractions]  # just differential fractions
  nTotFractions <- length(totProt)    # number of all fractions
  NSAfractions <- NSA[,1:nTotFractions]  # columns with NSA amounts; excludes additional columns


  # Compute amount of protein in a given fraction
  # yellow in Excel spreadsheet: "RSAtransformation instructions Lobel 29 Apr 2019"
  #amtProtFracOrig <- data.frame(markerProfiles %*% diag(totProt))  # old, inefficient way

  # use sweep operator
  # length(totProt)
  #  [1] 9
  # dim(markerProfiles)[1]
  #  [1] 8   # 8 rows
  # dim(markerProfiles)[2]
  # [1] 9    # 9 columns
  # Therefore, sweep across dimension 2

  AA <- data.frame(sweep(NSAfractions, 2,totProt, "*"))
  names(AA) <- colnames(NSAfractions)

  # Compute amount of given protein in fraction / amount of given protein in starting material
  # use these values to create mixtures
  # these are in red in Excel spreadsheet
  Acup <- data.frame(t(apply(AA, 1, function(x) x/sum(x[1:NstartMaterialFractions]))))
  names(Acup) <- colnames(AA)

  # relative specific activity (RSA)   # these are in heatmap in Excel spreadsheet
  #Difp <- sum(totProt[1:NstartMaterialFractions])   # total protein in the differential fractions
  #propFrac <- totProt/Difp  # proportion of protein in the differential fractions
  #rsa <- data.frame(as.matrix(relAmtProtFrac) %*% diag(1/propFrac))
  #names(rsa) <- colnames(markerProfiles)

  # return just Acup; no need for AA

 Acup
}



#' Compute relative specific amount from relative amounts in protein fractions
#' @param Acup amount of given protein in fraction / amount of that given protein in starting material
#' @param NstartMaterialFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param totProt Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's
#' @return rsa: relative specific amount
RSAfromAcup <- function(Acup, NstartMaterialFractions=6,
                         totProt=NULL) {

  if (ncol(Acup) != length(totProt)) {
    cat("Error from RSAfromAcup: no. of rows of Acup must match length of totProt\n")
  }
# # # # # # # #
# rename variabiles
# # # # # # # #

# Difp -> t.h
# propFrac -> t.cup

  # relative specific amount (RSA)   # these are in heatmap in Excel spreadsheet
  t.h <- sum(totProt[1:NstartMaterialFractions])   # total protein in the differential fractions
  t.cup <- totProt/t.h  # proportion of protein in the differential fractions (9 component vector)
  #rsa <- data.frame(as.matrix(relAmtProtFrac) %*% diag(1/propFrac))
  rsa <- sweep(Acup, 2, 1/t.cup, "*")
  names(rsa) <- colnames(Acup)
  rsa <- as.data.frame(rsa)
  # The following, which standardizes rsa rows to sum to one, returns the original markerProfiles !!
  NSAfromRSA <- t(apply(rsa,1, function(x) x/sum(x)))  # this is deprecated; no need

  rsa
}



#' compute a mixture of proteins in two compartments
#' @param AcupRef amount of given protein in fraction / amount of given protein in starting material
#' @param Loc1  row number of one compartment
#' @param Loc2  row number of other compartment
#' @param increment fraction increment from 0 to 1
#' @return mixAmount relative amounts of proteins in the fractions
proteinMix <- function(AcupRef, Loc1, Loc2, increment=0.10) {
  # replace mix.df with input.prop throughout
  # replace relAmount with Acup througout
  nrowRef <- nrow(AcupRef)
  if (Loc1 > Loc2) {
     Loc1orig <- Loc1
     Loc2orig <- Loc2
     Loc1 <- Loc2orig
     Loc2 <- Loc1orig
  }
  if (Loc2 > nrowRef) cat("Error, not enough rows\n")
  LocNames <- row.names(AcupRef)
  prop.vec <- seq(0,1,increment)
  qrop.vec <- 1 - prop.vec
  nrow.out <- length(prop.vec)
  Acup <- NULL
  mixProtNames <- NULL

  for (i in 1:nrow.out) {
    Acup.i <- prop.vec[i]*AcupRef[Loc1,] + qrop.vec[i]*AcupRef[Loc2,]
    Acup <- rbind(Acup, Acup.i)
    mixProtNames.i <- paste(prop.vec[i],"_", LocNames[Loc1], ":", qrop.vec[i], "_",LocNames[Loc2], sep='')
    mixProtNames <- c(mixProtNames, mixProtNames.i)

  }

  mixMat <- matrix(0, nrow=nrow.out, ncol=nrow(Acup))  # matrix of mixtures
  # each row is a "protein" with mixed residence
  # each column is a subcellular location, with the proportioal assignment
  mixMat[,Loc1] <- prop.vec
  mixMat[,Loc2] <- qrop.vec
  input.prop <- data.frame(mixMat)
  names(input.prop) <- row.names(Acup)
  rownames(input.prop) <- mixProtNames

  row.names(Acup) <- mixProtNames
  result <- list(Acup=Acup, input.prop=input.prop)
  result
  }


#' Directly compute relative specific activity from protProfileSummary (or from markerProfiles)
#'
#' @param protProfileLevels A matrix (with no prot names) giving the abundance level profiles of the subcellular locations, from 'cpaSetup', or many protein profiles
#' @param NstartMaterialFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param totProt Total protein counts in each of the fractions; this is necessary to compute RSA's
#'
#' @return protProfileRSA: relative specific activity

RSAfromNSA <- function(NSA=protProfileLevels, NstartMaterialFractions=6,
                      totProt=NULL) {

  missing.rows <- NSA[!complete.cases(NSA),]
  #if ({nrow(missing.rows) > 0} | {is.null(missing.rows)}) {
  #  cat("Error from RSAfromS: missing values not allowed\n")
  #  return(missing.rows)
  #}
  if (ncol(NSA) != length(totProt)) {
     cat("Error from rsaDirect: no. of rows of NSA must match length of totProt\n")
  }

  startMaterialFractions <- NSA[,1:NstartMaterialFractions]   # first columns of NSA
  nTotFractions <- length(totProt)    # number of all fractions
  NSAfractions <- NSA[,1:nTotFractions]  # columns with NSA amounts; excludes additional columns


  protAbund <- AcupFromNSA(NSAfractions,NstartMaterialFractions=NstartMaterialFractions,
                               totProt=totProt)
  Acup <- protAbund
  rsa <- RSAfromAcup(Acup, totProt=totProt)

  #Difp <- sum(totProt[1:NstartMaterialFractions])   # total protein in the differential fractions
  #TT <- as.matrix(protProfileLevels[,1:NstartMaterialFractions]) %*% totProt[1:NstartMaterialFractions]
  #protProfileRSA <- diag(as.numeric(1/TT)) %*% as.matrix(protProfileLevels) * Difp
  #protProfileRSA <- sweep(protProfileLevels, 1, 1/TT, "*")
  # truncate any large values at 15, for numerical stability

  rsa
 }


