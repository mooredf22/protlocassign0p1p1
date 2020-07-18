#' Set up reference profiles for constrained proportional assignment
#'
#' @param protProfileSummary data frame of protein identifiers (protName) and their relative abundance in centrifugation fractions.
#' @param refLocProteins List of reference proteins and their subcellular locations
#' @param n.channels Number of channels of abundance levels
#' @return A matrix markerLocR giving the abundance level profiles of the subcellular locations


cpaSetup <- function(protProfileSummary, refLocProteins=refLocProteins, n.channels) {
  # Find the mean profiles of the subcellular locations

  # "protProfileSummary" must be a data frame with these columns:
  #  Column 1: (protName): Prot name (or protein name)
  #  columns 2 - n.channels:   log2-transformed relative activities
  #  Last two columns:   Nspectra (number of spectra) and Nseq (number of distince peptide sequences)

  # "refLocProteins" must be a list of (1) the reference proteins and (2) the corresponding subcelluar locations

  # Select those entries from "protProfileSummary" that are reference proteins
  protNamesUpper <- toupper(protProfileSummary$protName)    # make names all upper case
  protProfileSummaryoriginal <- protProfileSummary
  protProfileSummary$protName <- protNamesUpper
  #meanCovarProtsStringent <- merge(x=refLocProteins, y=protProfileSummary,
  meanReferenceProts <- merge(x=refLocProteins, y=protProfileSummary,
         by.x="protName", by.y="protName", all.x=F, sort=F)

  # Find mean profiles for each sub-cellular location, using reference proteins
  #location.list <- unique(meanCovarProtsStringent$referenceCompartment)
  location.list <- as.character(unique(meanReferenceProts$referenceCompartment))
  n.loc <- length(location.list)
  meanProfile <- NULL
  for (i in 1:n.loc) {
    # i=1
    loc.i <- location.list[i]
    #profile.i <- meanCovarProtsStringent[meanCovarProtsStringent$referenceCompartment == loc.i,2+1:n.channels]
    profile.i <- meanReferenceProts[meanReferenceProts$referenceCompartment == loc.i,2+1:n.channels]
    meanProfile.i <- apply(profile.i, 2, mean, na.rm=T)
    meanProfile <- rbind(meanProfile, meanProfile.i)
    }
  row.names(meanProfile) <- location.list

  # Change name to "markerLoc" and its transpose "markerLocR" for compatibility with past
  # Then plot the profiles

  markerLoc <- t(meanProfile)
  markerLocR <- t(markerLoc)
  markerLocR
  }





Qfun4 <- function(pvec, y, gmat, methodQ="sumsquares") {
#   Assign sub-cellular location probabilities to
#   each protein. We use the "spg" function in package "BB".

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



protIndex <- function(protName, protProfileSummary=protProfileSummaryTMTms2,
                       exactMatch=F) {
  # return index of a protein name, or (if exactMatch=T) indices of proteins starting with
  #   the string given in "protName"
  n.prot <- nrow(protProfileSummary)

  prot.list <- toupper(as.character(protProfileSummary[,1]))  # must be in column 1
  #if (exactMatch) inx <- grep(paste("^", protName, "$", sep=""), prot.list, ignore.case=T)
  if (exactMatch) inx <- (1:n.prot)[protName == prot.list]
  if (!exactMatch) inx <- grep(paste("^",protName, sep=""), prot.list, ignore.case=T)
  if (length(inx) == 0) inx <- NA
  #if (protName %in% prot.list)  {inx <- (1:n.prot)[{prot.list == protName}]}
  #else if (protName %in% toupper(prot.list))  {inx <- (1:n.prot)[{toupper(prot.list) == protName}]}
  #is this all that is necessary?:
  #else if (toupper(protName) %in% toupper(prot.list))  {inx <- (1:n.prot)[{toupper(prot.list) == toupper(protName)}]}
  #else inx <- NA

  inx
  if (!is.na(inx[1])) {
    result <- data.frame(inx, prot.list[inx])
    names(result) <- c("Prot index number", "Prot name")
    }
  if (is.na(inx[1])) {
    result <- NA
    cat("protein not found\n")
  }
  result
  }
#protIndex("Tpp1", protProfileSummaryTMTms2)


protLocAssign <- function(i, protProfileSummary, markerLocR, n.channels, showProgress=T,
                          log2Transf=F, eps=2^(-5), maxit, assignProbsStart, maxR=NULL) {
  # maxit and assignPRobsStart must be specified
  # assignProbsStart must be NULL or have a column "protName" and assignment probabilities to use as starting values
  # use the spg function (in package BB) to assign proportionate assignments to compartments
  #nboot=1
  # i=2
  # i=1000
  if (!is.null(assignProbsStart)) {
    testEq <- {protProfileSummary$protName == assignProbsStart$protName}
    if (sum(as.numeric(testEq)) != nrow(protProfileSummary)) {
      cat("Error in protLocAssign: protName not identical in protProfileSummary and assignProbsStart\n")
    }
  }
  # check for excessively large component of markderLocR, if maxR is not NULL
  # If maxR is not NULL, truncate at maxR
  if (!is.null(maxR)) {
    markerLocRemp <- markerLocR
    ind.too.big <- {markerLocRtemp > maxR}
    markerLocR[ind.too.big] <- maxR
  }
  n.compartments <- nrow(markerLocR)
  allPeptideProfilesMat <- protProfileSummary[,1 + 1:n.channels]     # just the data

  yyT <- allPeptideProfilesMat[i,]



  yy <- as.numeric(yyT)       # leave them alone
  #yy <- 2^yyT - 0.04  # return to raw mode




  Nspectra.i <- protProfileSummary$Nspectra[i]   # this is the number of spectra for a protein
  Nseq.i <- protProfileSummary$Nseq[i] # number of unique sequences
  #channelsProb.i <- matrix(NA, nrow=nboot, ncol=8)

  startValsUnif <- rep(1/n.compartments, n.compartments)   # start with uniform probabilities
  if (!is.null(assignProbsStart)) startVals <- as.numeric(assignProbsStart[i, 2:(n.compartments + 1)]) # start with
  if (!log2Transf) {
   # first attempt at minimization: use uniform starting values
   temp <- try(BB::spg(startValsUnif, fn=Qfun4, project=proj.simplex, y=yy, gmat=t(markerLocR), methodQ="sumsquares", quiet=T,
    #temp <- try(spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(markerLocR), methodQ="sumabsvalue", quiet=T,
                                  control=list(maxit=maxit, trace=F)))
   convergeInd <- as.numeric((temp$message == "Successful convergence"))

   # If starting values are given (assignProbsStart is not null) and if the first attempt failed, try the starting values
   if ({convergeInd != 1} & {!is.null(assignProbsStart)}) {
   temp <- try(BB::spg(startVals, fn=Qfun4, project=proj.simplex, y=yy, gmat=t(markerLocR), methodQ="sumsquares", quiet=T,
                        #temp <- try(spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(markerLocR), methodQ="sumabsvalue", quiet=T,
                        control=list(maxit=maxit, trace=F)))
   }
  }
  if (log2Transf) {
    eps <- eps # defined in argument function
    markerLocRlog <- log2(markerLocR + eps)

    yyLog2= log2(yy + eps)
    temp <- try(BB::spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yyLog2, gmat=t(markerLocRlog), methodQ="sumsquares", quiet=T,
                    #temp <- try(spg(rep(1/n.compartments, n.compartments), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(markerLocR), methodQ="sumabsvalue", quiet=T,
                    control=list(maxit=maxit, trace=F)))
  }
  convergeInd <- as.numeric(!inherits(temp, "try-error"))

 # if (temp$message != "Successful convergence") {
 #   cat(paste(i, "\n"))
 #   cat(paste(as.character(yy),"\n"))
 #   cat(paste(as.character(protProfileSummary[i,])," \n"))
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
  if (convergeInd != 1) cat(paste("cpa does not converge for protein", i, "\n"))


  if (showProgress) {
    if (i == 500) print(paste(i, "proteins"))
    if ((i %% 1000) == 0) print(paste(i, "proteins"))
    }


    parEstTemp <- channelsMeanProb.i
    parEst <- c(parEstTemp, Nspectra.i, Nseq.i)


   parEst
   }

# test
#markerLocR <- cpaSetup(protProfileSummary, refLocProteins, n.channels=n.channels)
#protLocAssign(8, protProfileSummary, markerLocR, n.channels)
#protLocAssign(protIndex("Tpp1", protProfileSummary), protProfileSummary, markerLocR, n.channels)

# Test
if (F) {
yy= as.numeric(protProfileSummary[231,2:(n.channels+1)])
n.locs <- n.compartments
ans <- spg(rep(1/n.locs, n.locs), fn=Qfun4, project=proj.simplex, y=yy, gmat=t(markerLocR),
    control=list(trace=F))
 }

#' Carry out constrained proportional assignment on all proteinss
#'
#' @param protProfileSummary data frame of protein names and relative abundance levels.
#' @param markerLocR A matrix giving the abundance level profiles of the subcellular locations
#' @param n.channels Number of channels of abundance levels
#' @param log2Transf  Transform data before running CPA? Default is F
#' @param eps  Small quantity to add to data prior to taking a log transformation; ignored
#'              if log2Transf=F
#' @param maxit maximum number of iterations
#' @param assignProbsStart A matrix of starting values, one for each protein. The first column must be protName
#' @param maxR  A value to truncate values of markerLocR; if NULL (default), ignore
#' @return assignProbsOut  Data frame of proportionate assignments of each protein to compartments
proLocAll <- function(protProfileSummary, markerLocR, n.channels=n.channels, log2Transf=F, eps=2^(-5), maxit=10000,
                        assignProbsStart=NULL, maxR=NULL) {
  #markerLocR <- cpaSetup(protProfileSummary, refLocProteins, n.channels=n.channels)
  n.prot <- nrow(protProfileSummary)
  indList <- 1:n.prot
  #indList <- 1:50

   result <- sapply(indList, protLocAssign,
                protProfileSummary=protProfileSummary, markerLocR=markerLocR, n.channels=n.channels, showProgress=T,
                simplify=T, log2Transf=log2Transf, maxit=maxit, assignProbsStart=assignProbsStart,
                maxR=maxR)

  assignProbs <- data.frame(t(result))

  checkCols <- {ncol(assignProbs) == nrow(markerLocR) + 2}

  if (checkCols) {
    names(assignProbs) <- c(row.names(markerLocR), "Nspectra", "Npeptides")
    protName <- as.character(protProfileSummary$protName)  # make sure it is character
    assignProbsOut <- data.frame(protName, assignProbs)
    }
  if (!checkCols) {
    names(assignProbs) <- row.names(markerLocR)
    protName <- as.character(protProfileSummary$protName)
    assignProbsOut <- data.frame(protName, assignProbs)
  }
  assignProbsOut
  }


proLocAllmcore <- function(protProfileSummary, refLocProteins, n.channels=n.channels) {
  # experimental multicore version
  markerLocR <- cpaSetup(protProfileSummary, refLocProteins, n.channels=n.channels)
  n.prot <- nrow(protProfileSummary)
  indList <- 1:n.prot
#  result <- sapply(indList, protLocAssign,
#                protProfileSummary=protProfileSummary, markerLocR=markerLocR, n.channels=7, showProgress=T,
#                simplify=T)
  result <- sfLapply(indList, protLocAssign,
                protProfileSummary=protProfileSummary, markerLocR=markerLocR, n.channels=7)
  assignProbs <- data.frame(t(result))
  #dim(assignProbs)
  # [1] 278867      8
  #dim(allPeptideProfilesMat)
  # [1] 278867      7
  names(assignProbs) <- c(row.names(markerLocR), "Nspectra", "Npeptides")
  protName <- as.characther(protProfileSummary$protName)
  assignProbsOut <- data.frame(protName, assignProbs)
  assignProbsOut
  }


protCap <- function(x) {
    # convert single protein name to capitalize first character only
    # Typically, human proteins are all in capitals and
    #   rat are capitalized (first character only) with the rest in lower case.
    firstLetter <- toupper(substring(x, 1, 1))
    remainingLetters <- tolower(substring(x,2))
    paste(firstLetter, remainingLetters, sep='')
}

#' assign proteins to a subcellular location using CPA estimates
#'
#' @param assignLocProps matrix of proportion estimates for each protein
#' @param cutoff cutoff for assigning a protein to a location
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
