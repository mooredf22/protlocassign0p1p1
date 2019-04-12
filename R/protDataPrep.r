#' Prepare data for constrained proportional assignment by reweighting differential and nyc fractions
#'
#' @param geneProfileSummary data frame of gene/protein identifiers (geneName) and their relative abundance in centrifugation fractions.
#' @param n.diff.chan Number of channels of differential centrifugation fractions
#' @param n.nyc.chan Number of nycodenz fractions (derived from the L1 fraction)
#' @return A matrix geneProfileSummary with re-weighted relative abundance fractions

#'

protDataPrep <- function(geneProfileSummary, n.diff.chan=6, n.nyc.chan=3) {
   # prepare the geneProfileSummary data with different constraints
   # for differential fractions and Nyc fractions

  # "geneProfileSummary" must be a data frame with these columns:
  #  Column 1: (geneName): Gene name (or protein name)
  #  columns 2 - n.channels:   log2-transformed relative activities
  #  Last two columns:   Nspectra (number of spectra) and Nseq (number of distince peptide sequences)

  # The differential fractions are N, M, L1, L2, P, S   (n.diff.chan = 6, ususally)
  # L1 must be the third fraction
  # The Nyc fractions are Nyc.1, Nyc.2, Nyc.3, or just Nyc.2, or absent

  # Note: the Nyc fractions are derived from the L1 fraction

  # options:
  #   n.nyc.chan = 3:  (1) constrain the differential fractions to sum to one
  #                   (2) constrain the Nyc fractions to sum to one
  #                   (3) multiply the Nyc fractions by the L1 amount

  #   n.nyc.chan = 1   (1) constrain all fractions to sum to one
  #                   (2) multiply the Nyc.2 fraction by L1
  #


  if (n.nyc.chan >= 2) {
    allFractionsT <- geneProfileSummary[, 2:(1 + n.diff.chan + n.nyc.chan)]
    differentialFractionsT <- allFractionsT[,1:n.diff.chan]
    nycFractionsT <- allFractionsT[,(1 + n.diff.chan):(n.diff.chan + n.nyc.chan)]
    differentialFractions <- t(apply(differentialFractionsT,1, function(x) x/sum(x)))  # constrain rows to sum to 1
    nycFractions <- t(apply(nycFractionsT, 1, function(x) x/sum(x)))*differentialFractions[,3]   # constrain rows to sum to 1; mult by L1
    allFractionsT <- data.frame(differentialFractions, nycFractions)
  }
  if (n.nyc.chan == 1) {
    allFractionsT <- geneProfileSummary[, 2:(1 + n.diff.chan + n.nyc.chan)]
    differentialFractionsT <- allFractionsT[,1:n.diff.chan]
    differentialFractions <- t(apply(differentialFractionsT,1, function(x) x/sum(x)))  # constrain rows to sum to 1
    nycFractions <- allFractionsT[,(1 + n.diff.chan):(n.diff.chan + n.nyc.chan)]*differentialFractions[,3]  #mult by normalized L1
    allFractionsT <- data.frame(differentialFractions, nycFractions)
  }
  if (n.nyc.chan == 0) {
    allFractionsT <- geneProfileSummary[, 2:(1 + n.diff.chan + n.nyc.chan)]
    allFractions <- t(apply(allFractionsT, 1, function(x) x/sum(x)))

  }
  allFractions <- t(apply(allFractionsT, 1, function(x) x/sum(x)))  # constrain to sum to one, after all previous transformations

  geneProfileSummary[,2:(1 + n.diff.chan + n.nyc.chan)] <- allFractions
  geneProfileSummary

}



