

#' Set up RSA (relative specific amount) profiles for constrained proportional assignment
#'
#' @param SS A matrix giving the specific amount or normalized specific amount
#'           profiles, either marker profiles from 'cpaSetup' (which requires
#'           normalized specific amounts as input), or a list of many protein profiles
#' @param nDiffFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param nNycFractions Number of Nycodenz fractions, typically 3, but could be 1 (if only Nyc2 is present) or 0 if none
#' @param totProt Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's
#'
#' @return amtProtFrac: amount of protein in a given fraction
#' @return relAmtProtFrac amount of given protein in fraction / amount of given protein in starting material

# # # # # # # # #
# name changes July 8, 2020
# # # # # # # # #

# function name change
#abundanceTransform <- function(markerLocR, nDiffFractions=6, nNycFractions=3,

# markerLocR -> SS   # no need to limit this to the marker profiles
# amtProtFrac -> AA  # these are the amounts of protein (or other species)
# relAmtProtFrac -> Acup  # relative amounts
                # normalize specific amounts in each fraction
                # to the amount in starting material
    # these are of interest for simulations
#

relAmtTransform <- function(SS, nDiffFractions=6, nNycFractions=3,
                               totProt=NULL) {
  if ((nDiffFractions + nNycFractions) != ncol(SS)) {
    cat("Error from RSAtransform\nTotal number of fractions must be the number of columns of SS\n")
    return()
  }
  diffFractions <- SS[,1:nDiffFractions]
  if (nNycFractions > 0) nycFractions <- SS[,(nDiffFractions+1):(nDiffFractions+nNycFractions)]
  if (nNycFractions == 0) nycFractions <- NULL

  # Compute amount of protein in a given fraction
  # yellow in Excel spreadsheet: "RSAtransformation instructions Lobel 29 Apr 2019"
  #amtProtFracOrig <- data.frame(markerLocR %*% diag(totProt))  # old, inefficient way

  # use sweep operator
  # length(totProt)
  #  [1] 9
  # dim(markerLocR)[1]
  #  [1] 8   # 8 rows
  # dim(markerLocR)[2]
  # [1] 9    # 9 columns
  # Therefore, sweep across dimension 2

  AA <- data.frame(sweep(SS, 2,totProt, "*"))
  names(AA) <- colnames(SS)

  # Compute amount of given protein in fraction / amount of given protein in starting material
  # use these values to create mixtures
  # these are in red in Excel spreadsheet
  Acup <- data.frame(t(apply(AA, 1, function(x) x/sum(x[1:nDiffFractions]))))
  names(Acup) <- colnames(AA)

  # relative specific activity (RSA)   # these are in heatmap in Excel spreadsheet
  #Difp <- sum(totProt[1:nDiffFractions])   # total protein in the differential fractions
  #propFrac <- totProt/Difp  # proportion of protein in the differential fractions
  #rsa <- data.frame(as.matrix(relAmtProtFrac) %*% diag(1/propFrac))
  #names(rsa) <- colnames(markerLocR)

  # The following, which standardizes rsa rows to sum to one, returns the original markerLocR !!
  #rsaFractions <- t(apply(rsa,1, function(x) x/sum(x)))

  result <- list(AA=AA, Acup=Acup)
  result
}



#' Compute relative specific activity and RSA fractions from rlatve amounts in protein fractions
#' @param Acup amount of given protein in fraction / amount of that given protein in starting material
#' @param nDiffFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param nNycFractions Number of Nycodenz fractions, typically 3, but could be 1 (if only Nyc2 is present) or 0 if none
#' @param totProt Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's
#' @return rsa: relative specific amount
RSAtransform <- function(Acup, nDiffFractions=6, nNycFractions=3,
                         totProt=NULL) {
# # # # # # # #
# rename variabiles
# # # # # # # #

# Difp -> t.h
# propFrac -> t.cup

  # relative specific amount (RSA)   # these are in heatmap in Excel spreadsheet
  t.h <- sum(totProt[1:nDiffFractions])   # total protein in the differential fractions
  t.cup <- totProt/t.h  # proportion of protein in the differential fractions (9 component vector)
  #rsa <- data.frame(as.matrix(relAmtProtFrac) %*% diag(1/propFrac))
  rsa <- sweep(Acup, 2, 1/t.cup, "*")
  names(rsa) <- colnames(Acup)

  # The following, which standardizes rsa rows to sum to one, returns the original markerLocR !!
  rsaFractions <- t(apply(rsa,1, function(x) x/sum(x)))  # this is deprecated; no need

  rsa
}

#RSAtransform(relAmtProtFrac)

#RSAtransform(mixCytoLyso)

#' compute a mixture of proteins in two compartments
#' @param relAmtProtFrac amount of given protein in fraction / amount of given protein in starting material
#' @param Loc1  row number of one compartment
#' @param Loc2  row number of other compartment
#' @param increment fraction increment from 0 to 1
#' @return mixAmount relative amounts of proteins in the fractions
proteinMix <- function(Acup, Loc1, Loc2, increment=0.10) {
  nrowRef <- nrow(Acup)
  if (Loc1 > Loc2) {
     Loc1orig <- Loc1
     Loc2orig <- Loc2
     Loc1 <- Loc2orig
     Loc2 <- Loc1orig
  }
  if (Loc2 > nrowRef) cat("Error, not enough rows\n")
  LocNames <- row.names(Acup)
  prop.vec <- seq(0,1,increment)
  qrop.vec <- 1 - prop.vec
  nrow.out <- length(prop.vec)
  mixAmount <- NULL
  mixProtNames <- NULL
  for (i in 1:nrow.out) {
    mixAmount.i <- prop.vec[i]*Acup[Loc1,] + qrop.vec[i]*Acup[Loc2,]
    mixAmount <- rbind(mixAmount, mixAmount.i)
    mixProtNames.i <- paste(prop.vec[i],"_", LocNames[Loc1], ":", qrop.vec[i], "_",LocNames[Loc2], sep='')
    mixProtNames <- c(mixProtNames, mixProtNames.i)
  }
  row.names(mixAmount) <- mixProtNames
  mixAmount
  }

#mixCytoLyso <- proteinMix(relAmtProtFrac, 1, 4)
#RSAtransform(mixCytoLyso)$rsaFractions

#' Directly compute relative specific activity from protProfileSummary (or from markerLocR)
#'
#' @param protProfileLevels A matrix (with no prot names) giving the abundance level profiles of the subcellular locations, from 'cpaSetup', or many protein profiles
#' @param nDiffFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param nNycFractions Number of Nycodenz fractions, typically 3, but could be 1 (if only Nyc2 is present) or 0 if none
#' @param totProt Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's
#'
#' @return protProfileRSA: relative specific activity

RSAfromS <- function(SS=protProfileLevels, nDiffFractions=6, nNycFractions=3,
                      totProt=NULL) {

  missing.rows <- SS[!complete.cases(SS),]
  if ({nrow(missing.rows) > 0} | {is.null(missing.rows)}) {
    cat("Error from rsaDirect: missing values not allowed\n")
    return(missing.rows)
  }
  if ((nDiffFractions + nNycFractions) != ncol(SS)) {
    cat("Error from RSAdirect\nTotal number of fractions must be the number of columns of markerLocR\n")
    return()
  }
  diffFractions <- SS[,1:nDiffFractions]
  if (nNycFractions > 0) nycFractions <- SS[,(nDiffFractions+1):(nDiffFractions+nNycFractions)]
  if (nNycFractions == 0) nycFractions <- NULL

  protAbund <- relAmtTransform(SS,nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
  Acup <- protAbund$Acup
  rsa <- RSAtransform(Acup, totProt=totProtUse)

  #Difp <- sum(totProt[1:nDiffFractions])   # total protein in the differential fractions
  #TT <- as.matrix(protProfileLevels[,1:nDiffFractions]) %*% totProt[1:nDiffFractions]
  #protProfileRSA <- diag(as.numeric(1/TT)) %*% as.matrix(protProfileLevels) * Difp
  #protProfileRSA <- sweep(protProfileLevels, 1, 1/TT, "*")
  # truncate any large values at 15, for numerical stability

  rsa
 }

# markerLocRrsa <- rsaDirect(protProfileLevels=markerLocR, totProt=c(46.044776, 48.955954, 1.384083, 1.566324, 24.045584, 58.181818, 0.0684596))
