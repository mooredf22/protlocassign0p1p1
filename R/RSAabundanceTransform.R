

#' Set up RSA (relative specific activity) profiles for constrained proportional assignment
#'
#' @param matLocR A matrix giving the abundance level profiles of the subcellular locations, from 'cpaSetup', or many gene profiles
#' @param nDiffFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param nNycFractions Number of Nycodenz fractions, typically 3, but could be 1 (if only Nyc2 is present) or 0 if none
#' @param totProt Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's
#'
#' @return amtProtFrac: amount of protein in a given fraction

abundanceTransform <- function(matLocR, nDiffFractions=6, nNycFractions=3,
                               totProt=c(46.044776, 48.955954, 1.384083, 1.566324, 24.045584, 58.181818, 0.0368564, 0.0684596, 1.27301587)) {
  if ((nDiffFractions + nNycFractions) != ncol(matLocR)) {
    cat("Error from RSAtransform\nTotal number of fractions must be the number of columns of matLocR\n")
    return()
  }
  diffFractions <- matLocR[,1:nDiffFractions]
  if (nNycFractions > 0) nycFractions <- matLocR[,(nDiffFractions+1):(nDiffFractions+nNycFractions)]
  if (nNycFractions == 0) nycFractions <- NULL

  # Compute amount of protein in a given fraction   # yellow in Excel spreadsheet: "RSAtransformation instructions Lobel 29 Apr 2019"
  amtProtFrac <- matLocR %*% diag(totProt)

  # Compute amount of given protein in fraction / amount of given protein in starting material
  # use these values to create mixtures
  # these are in red in Excel spreadsheet
  relAmtProtFrac <- data.frame(t(apply(amtProtFrac, 1, function(x) x/sum(x[1:nDiffFractions]))))
  names(relAmtProtFrac) <- colnames(matLocR)

  # relative specific activity (RSA)   # these are in heatmap in Excel spreadsheet
  #Difp <- sum(totProt[1:nDiffFractions])   # total protein in the differential fractions
  #propFrac <- totProt/Difp  # proportion of protein in the differential fractions
  #rsa <- data.frame(as.matrix(relAmtProtFrac) %*% diag(1/propFrac))
  #names(rsa) <- colnames(matLocR)

  # The following, which standardizes rsa rows to sum to one, returns the original matLocR !!
  #rsaFractions <- t(apply(rsa,1, function(x) x/sum(x)))

  result <- list(amtProtFrac=amtProtFrac, relAmtProtFrac=relAmtProtFrac)
  result
}

#tempRel <- abundanceTransform(matLocR)
#relAmtProtFrac <- tempRel$relAmtProtFrac

#mixCytoLyso10 <- 0.1*relAmtProtFrac[1,] + 0.9*relAmtProtFrac[4,]
#mixCytoLyso50 <- 0.5*relAmtProtFrac[1,] + 0.5*relAmtProtFrac[4,]
#mixCytoLyso90 <- 0.9*relAmtProtFrac[1,] + 0.1*relAmtProtFrac[4,]
#mixCytoLyso <- rbind(mixCytoLyso10, mixCytoLyso50, mixCytoLyso90)

#' @param relAmtProtFrac amount of given protein in fraction / amount of given protein in starting material
#' @param nDiffFractions Number of differential fractions, typically 6, for N, M, L1, L2, P, and S
#' @param nNycFractions Number of Nycodenz fractions, typically 3, but could be 1 (if only Nyc2 is present) or 0 if none
#' @param totProt Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's
#' @return rsa: relative specific activity
RSAtransform <- function(relAmtProtFrac, nDiffFractions=6, nNycFractions=3,
                         totProt=c(46.044776, 48.955954, 1.384083, 1.566324, 24.045584, 58.181818, 0.0368564, 0.0684596, 1.27301587)) {

  # relative specific activity (RSA)   # these are in heatmap in Excel spreadsheet
  Difp <- sum(totProt[1:nDiffFractions])   # total protein in the differential fractions
  propFrac <- totProt/Difp  # proportion of protein in the differential fractions
  rsa <- data.frame(as.matrix(relAmtProtFrac) %*% diag(1/propFrac))
  names(rsa) <- colnames(matLocR)

  # The following, which standardizes rsa rows to sum to one, returns the original matLocR !!
  rsaFractions <- t(apply(rsa,1, function(x) x/sum(x)))

  result <- list(rsa=rsa, rsaFractions=rsaFractions)
  result
}

#RSAtransform(relAmtProtFrac)

#RSAtransform(mixCytoLyso)


