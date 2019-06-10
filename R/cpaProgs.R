#' Set up reference profiles for constrained proportional assignment
#'
#' @param geneProfileSummary data frame of gene/protein identifiers (geneName) and their relative abundance in centrifugation fractions.
#' @param refLocProteins List of reference proteins and their subcellular locations
#' @param n.channels Number of channels of abundance levels
#' @return A matrix matLocR giving the abundance level profiles of the subcellular locations

#'

cpaSetup <- function(geneProfileSummary, refLocProteins=refLocProteins, n.channels) {
  # Find the mean profiles of the subcellular locations

  # "geneProfileSummary" must be a data frame with these columns:
  #  Column 1: (geneName): Gene name (or protein name)
  #  columns 2 - n.channels:   log2-transformed relative activities
  #  Last two columns:   Nspectra (number of spectra) and Nseq (number of distince peptide sequences)

  # "refLocProteins" must be a list of (1) the reference proteins and (2) the corresponding subcelluar locations

  # Select those entries from "geneProfileSummary" that are reference genes
  geneNamesUpper <- toupper(geneProfileSummary$geneName)    # make names all upper case
  geneProfileSummaryoriginal <- geneProfileSummary
  geneProfileSummary$geneName <- geneNamesUpper
  #meanCovarGenesStringent <- merge(x=refLocProteins, y=geneProfileSummary,
  meanReferenceGenes <- merge(x=refLocProteins, y=geneProfileSummary,
         by.x="geneName", by.y="geneName", all.x=F, sort=F)

  # Find mean profiles for each sub-cellular location, using reference genes
  #location.list <- unique(meanCovarGenesStringent$referenceCompartment)
  location.list <- unique(meanReferenceGenes$referenceCompartment)
  n.loc <- length(location.list)
  meanProfile <- NULL
  for (i in 1:n.loc) {
    # i=1
    loc.i <- location.list[i]
    #profile.i <- meanCovarGenesStringent[meanCovarGenesStringent$referenceCompartment == loc.i,2+1:n.channels]
    profile.i <- meanReferenceGenes[meanReferenceGenes$referenceCompartment == loc.i,2+1:n.channels]
    meanProfile.i <- apply(profile.i, 2, mean, na.rm=T)
    meanProfile <- rbind(meanProfile, meanProfile.i)
    }
  row.names(meanProfile) <- location.list

  # Change name to "matLoc" and its transpose "matLocR" for compatibility with past
  # Then plot the profiles

  matLoc <- t(meanProfile)
  matLocR <- t(matLoc)
  matLocR
  }





Qfun4 <- function(pvec, y, gmat, methodQ="sumsquares") {
#   Assign sub-cellular location probabilities to
#   each gene. We use the "spg" function in package "BB".

   resultA <- y - pvec%*%t(gmat)

   if (methodQ=="sumsquares") result <- sum(resultA^2)
   if (methodQ=="sumabsvalue") result <- sum(abs(resultA))
   result
   }

#library(BB)
proj.simplex <- function(y, c=1) {
  #####
  # project an n-dim vector y to the simplex S_n
  # S_n = { x | x \in R^n, 0 <= x <= c, sum(x) = c}
  # Copyright:  Ravi Varadhan, Johns Hopkins University
  # August 8, 2012

  #####

  n <- length(y)
  sy <- sort(y, decreasing=TRUE)
  csy <- cumsum(sy)
  rho <- max(which(sy > (csy - c)/(1:n)))
  theta <- (csy[rho] - c) / rho
  return(pmax(0, y - theta))
  }



protIndex <- function(protName, geneProfileSummary) {
  # return index of a protein name
  n.prot <- nrow(geneProfileSummary)

  prot.list <- as.character(geneProfileSummary[,1])  # must be in column 1
  if (protName %in% prot.list)  {inx <- (1:n.prot)[{prot.list == protName}]}
  else if (protName %in% toupper(prot.list))  {inx <- (1:n.prot)[{toupper(prot.list) == protName}]}
  #is this all that is necessary?:
  else if (toupper(protName) %in% toupper(prot.list))  {inx <- (1:n.prot)[{toupper(prot.list) == toupper(protName)}]}
  else inx <- NA

  inx
  }
#protIndex("Tpp1", geneProfileSummary)
#protIndex("TPP1", geneProfileSummary)


protLocAssign <- function(i, geneProfileSummary, matLocR, n.channels, showProgress=T, log2Transf=F, maxit, assignProbsStart) {
  # maxit and assignPRobsStart must be specified
  # assignProbsStart must be NULL or have a column "geneName" and assignment probabilities to use as starting values
  # use the spg function (in package BB) to assign proportionate assignments to compartments
  #nboot=1
  # i=2
  # i=1000
  if (!is.null(assignProbsStart)) {
    testEq <- {geneProfileSummary$geneName == assignProbsStart$geneName}
    if (sum(as.numeric(testEq)) != nrow(geneProfileSummary)) {
      cat("Error in protLocAssign: geneName not identical in geneProfileSummary and assignProbsStart\n")
    }
  }
  n.compartments <- nrow(matLocR)
  allPeptideProfilesMat <- geneProfileSummary[,1 + 1:n.channels]     # just the data

  yyT <- allPeptideProfilesMat[i,]



  yy <- as.numeric(yyT)       # leave them alone
  #yy <- 2^yyT - 0.04  # return to raw mode




  Nspectra.i <- geneProfileSummary$Nspectra[i]   # this is the number of spectra for a gene
  Nseq.i <- geneProfileSummary$Nseq[i] # number of unique sequences
  #channelsProb.i <- matrix(NA, nrow=nboot, ncol=8)

  startValsUnif <- rep(1/n.compartments, n.compartments)   # start with uniform probabilities
  if (!is.null(assignProbsStart)) startVals <- as.numeric(assignProbsStart[i, 2:(n.compartments + 1)]) # start with
  if (!log2Transf) {
   # first attempt at minimization: use uniform starting values
   temp <- try(BB::spg(startValsUnif, fn=Qfun4, project=proj.simplex, y=yy, gmat=t(matLocR), methodQ="sumsquares", quiet=T,
    #temp <- try(spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(matLocR), methodQ="sumabsvalue", quiet=T,
                                  control=list(maxit=maxit, trace=F)))
   convergeInd <- as.numeric((temp$message == "Successful convergence"))

   # If starting values are given (assignProbsStart is not null) and if the first attempt failed, try the starting values
   if ({convergeInd != 1} & {!is.null(assignProbsStart)}) {
   temp <- try(BB::spg(startVals, fn=Qfun4, project=proj.simplex, y=yy, gmat=t(matLocR), methodQ="sumsquares", quiet=T,
                        #temp <- try(spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(matLocR), methodQ="sumabsvalue", quiet=T,
                        control=list(maxit=maxit, trace=F)))
   }
  }
  if (log2Transf) {
    eps <- 2^(-5)
    matLocRlog <- log2(matLocR + eps)

    yyLog2= log2(yy + eps)
    temp <- try(BB::spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yyLog2, gmat=t(matLocRlog), methodQ="sumsquares", quiet=T,
                    #temp <- try(spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(matLocR), methodQ="sumabsvalue", quiet=T,
                    control=list(maxit=maxit, trace=F)))
  }
  convergeInd <- as.numeric(!inherits(temp, "try-error"))
 # if (temp$message != "Successful convergence") {
 #   cat(paste(i, "\n"))
 #   cat(paste(as.character(yy),"\n"))
 #   cat(paste(as.character(geneProfileSummary[i,])," \n"))
 # }
  channelsMeanProb.i <- rep(NA, n.compartments)
  #channelsProbCImat.i <- matrix(NA, nrow=2, ncol=8)
  convCrit.i <- NA
  feval.i <- NA
  if (!is.atomic(temp)) {
     channelsMeanProb.i <- temp$par
     convCrit.i <- temp$gradient
     feval.i <- temp$feval
     convergeInd <- as.numeric((temp$message == "Successful convergence"))
     }
  nNoConverge.i <- 0



  if (showProgress) {
    if ((i %% 500) == 0) print(paste(i, "genes"))
    }


    parEstTemp <- channelsMeanProb.i
    parEst <- c(parEstTemp, Nspectra.i, Nseq.i, convergeInd)


   parEst
   }

# test
#matLocR <- cpaSetup(geneProfileSummary, refLocProteins, n.channels=n.channels)
#protLocAssign(8, geneProfileSummary, matLocR, n.channels)
#protLocAssign(protIndex("Tpp1", geneProfileSummary), geneProfileSummary, matLocR, n.channels)

# Test
if (F) {
yy= as.numeric(geneProfileSummary[231,2:(n.channels+1)])
n.locs <- n.compartments
ans <- spg(rep(1/n.locs, n.locs), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(matLocR),
    control=list(trace=F))
 }

#' Carry out constrained proportional assignment on all genes
#'
#' @param geneProfileSummary data frame of protein names and relative abundance levels.
#' @param matLocR A matrix giving the abundance level profiles of the subcellular locations
#' @param n.channels Number of channels of abundance levels
#' @param maxit maximum number of iterations
#' @param assignProbsStart A matrix of starting values, one for each gene. The first column must be geneName
#' @return assignProbsOut  Data frame of proportionate assignments of each protein to compartments
proLocAll <- function(geneProfileSummary, matLocR, n.channels=n.channels, log2Transf=F, maxit=10000, assignProbsStart=NULL) {
  #matLocR <- cpaSetup(geneProfileSummary, refLocProteins, n.channels=n.channels)
  n.gene <- nrow(geneProfileSummary)
  indList <- 1:n.gene
  #indList <- 1:50

   result <- sapply(indList, protLocAssign,
                geneProfileSummary=geneProfileSummary, matLocR=matLocR, n.channels=n.channels, showProgress=T,
                simplify=T, log2Transf=log2Transf, maxit=maxit, assignProbsStart=assignProbsStart)

  assignProbs <- data.frame(t(result))
  #dim(assignProbs)
  # [1] 278867      8
  #dim(allPeptideProfilesMat)
  # [1] 278867      7
  checkCols <- {ncol(assignProbs) == nrow(matLocR) + 3}
  if (checkCols) {
    names(assignProbs) <- c(row.names(matLocR), "Nspectra", "Npeptides", "convergeInd")
    geneName <- geneProfileSummary$geneName
    assignProbsOut <- data.frame(geneName, assignProbs)
    }
  if (!checkCols) {
    names(assignProbs) <- c(row.names(matLocR), "convergeInd")
    geneName <- geneProfileSummary$geneName
    assignProbsOut <- data.frame(geneName, assignProbs)
  }
  assignProbsOut
  }


proLocAllmcore <- function(geneProfileSummary, refLocProteins, n.channels=n.channels) {
  # experimental multicore version
  matLocR <- cpaSetup(geneProfileSummary, refLocProteins, n.channels=n.channels)
  n.gene <- nrow(geneProfileSummary)
  indList <- 1:n.gene
#  result <- sapply(indList, protLocAssign,
#                geneProfileSummary=geneProfileSummary, matLocR=matLocR, n.channels=7, showProgress=T,
#                simplify=T)
  result <- sfLapply(indList, protLocAssign,
                geneProfileSummary=geneProfileSummary, matLocR=matLocR, n.channels=7)
  assignProbs <- data.frame(t(result))
  #dim(assignProbs)
  # [1] 278867      8
  #dim(allPeptideProfilesMat)
  # [1] 278867      7
  names(assignProbs) <- c(row.names(matLocR), "Nspectra", "Npeptides")
  geneName <- geneProfileSummary$geneName
  assignProbsOut <- data.frame(geneName, assignProbs)
  assignProbsOut
  }


geneCap <- function(x) {
    # convert single gene name to capitalize first character only
    # Typically, human genes are all in capitals and
    #   rat are capitalized (first character only) with the rest in lower case.
    firstLetter <- toupper(substring(x, 1, 1))
    remainingLetters <- tolower(substring(x,2))
    paste(firstLetter, remainingLetters, sep='')
}

#' assign genes to a subcellular location using CPA estimates
#'
#' @param assignLocProps matrix of proportion estimates for each gene
#' @param cutoff cutoff for assigning a gene to a location
#' @param Locations list of subcellular locations

assignCPAloc <- function(assignLocProps, cutoff=0.8, Locations) {
  # assign location to the one with a proportion > cutoff
  # assignLocProps <- assignPropsAll[1,]
  if (cutoff <= 0.5) cutoff <- 0.5001
  exceedCutoff <- {assignLocProps > cutoff}
  if (sum(exceedCutoff) > 0) {
    indLoc <- which(exceedCutoff)
    protLoc <- Locations[indLoc]
  }
  if (sum(exceedCutoff) == 0) protLoc <- "unclassified"
  protLoc
}
