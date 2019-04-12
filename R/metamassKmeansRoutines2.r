#' Assign genes to locations using metamass procedure
#'
#' @param geneProfileSummary data frame of protein names and relative abundance levels.
#' @param refLocProteins List of reference proteins and their subcellular locations
#' @param nCol  number of columns of abundance levels
#' @param clustSize  average cluster size
#' @return assignProbsOut  assignments of each protein to compartments

protCluster <- function(geneProfileSummary, refLocProteins, refCol=1, nCol=9, clustSize=5) {
  # geneProfileSummary must contain gene names in first column, abundance levels
  #   in next nCol columns, all positive and summing to 1,
  #   nCol is the number of channels
  #   clustSize is the average cluster size
  #
  nColTot <- refCol + nCol
  lookup <- refLocProteins
  names(lookup) <- c("gene", "compartment")

  geneMatch <- function(x) as.character(lookup$compartment[match(x, lookup$gene)])
  #geneMatch("AARS")
  #geneMatch("AARS2")

  classif <- unlist(lapply(geneProfileSummary$geneName, geneMatch))
  table(classif)

  newProtClass <- data.frame(geneProfileSummary, classif) # this data frame has the "overlap" classification
  dim(newProtClass)
  # [1] 857064     14
  #refCol <- refCol
  exptData <- newProtClass[,(refCol + 1):(refCol + nCol)]        # just the data
  #clustSize <- 5
  #clustSize <- 10
  n.clusters <- round(nrow(exptData)/clustSize)
  n.clusters        # average of 5 per cluster
  # [1] 789
  result.kmeans <- kmeans(as.matrix(exptData), centers=n.clusters)
  #plot(result.kmeans)  # doesn't work

  cl <- result.kmeans$cluster
  newProtClassCl <- data.frame(cl, newProtClass)
#data.frame(1:length(names(newProtClass)),names(newProtClass))
#clProts <- data.frame(cl, newProtClass[,c(1, 2, 3:12)])
  clOrdInd <- order(newProtClassCl$cl)
  clProtsOrd <- newProtClassCl[clOrdInd, ]        # ordered by cluster number

  #table(clProtsOrd$classif)
  # Now look through all of the clusters. For each, find the
  #   most common location, and assign all to that location.

  purity.vec <- NULL          # number of genes classifed correctly / number of genes classified
  maxCat.vec <- NULL          # most common classification
  n.clust.vec <- NULL         # cluster size
  num.classified.vec <- NULL  # number of genes classified
  purity.mat <- NULL          # matrix of purities by classification category
  for (i in 1:n.clusters) {
    # i=1
    # i=2
    clProtsOrd.i <- clProtsOrd[clProtsOrd$cl == i,]
    temp <- table(clProtsOrd.i$classif)
    tableValues <- as.numeric(temp)
    tableNames <- names(temp)         # names of subcellular locations
    wh.max <- 0
    num.classified.i <- 0     # number of genes in a cluster with classifications; start with 0
    if (sum(tableValues) > 0) {
      wh.max <- which.max(tableValues)
      wh.max.rev <- which.max(rev(tableValues))
      if (wh.max + wh.max.rev == (length(tableValues) + 1))  {       # true if the max value is unique
        num.classified.i <- sum(tableValues)   #number of genes in cluster that have been classified
        purity.i <- tableValues[wh.max]/num.classified.i   # purity: proportion classified as maximum
        purity.vec <- c(purity.vec, rep(purity.i, nrow(clProtsOrd.i)))
        maxCat <- tableNames[wh.max]
        }
      if (wh.max + wh.max.rev != (length(tableValues) + 1))  {   #max value not unique
        purity.i <- NA
        purity.vec <- c(purity.vec, rep(purity.i, nrow(clProtsOrd.i)))
        maxCat <- "none"
        }
      }
    if (sum(tableValues) == 0) {
       maxCat <- "unclassified"
       purity.i <- NA
       purity.vec <- c(purity.vec, rep(purity.i, nrow(clProtsOrd.i)))
       }
    maxCat.vec <- c(maxCat.vec, rep(maxCat, nrow(clProtsOrd.i)))
    num.classified.vec <- c(num.classified.vec, rep(num.classified.i, nrow(clProtsOrd.i)))
    n.clust.i <- nrow(clProtsOrd.i)
    n.clust.vec <- c(n.clust.vec, rep(n.clust.i, nrow(clProtsOrd.i)))
    for (jj in 1:n.clust.i) {
      if (num.classified.i > 0) purity.mat <- rbind(purity.mat, tableValues/num.classified.i)
      if (num.classified.i == 0) purity.mat <- rbind(purity.mat, tableValues)
    }
  }
  purity.df <- data.frame(purity.mat)
  names(purity.df) <- tableNames
  geneName <- clProtsOrd$geneName

  purityStats <- data.frame(purity.vec, maxCat.vec, num.classified.vec, n.clust.vec)

  row.names(clProtsOrd) <- geneName
  row.names(purityStats) <- geneName
  row.names(purity.df) <- geneName

   clProtsOrdAll <- list(clProtsOrd=clProtsOrd, purityStats=purityStats, purity.df=purity.df)
   clProtsOrdAll
  }


#' Assign genes to locations using metamass procedure; Carries out protCluster multiple
#'   times to assess consistency of assignments to compartments
#'
#' @param geneProfileSummary data frame of protein names and relative abundance levels.
#' @param refLocProteins List of reference proteins and their subcellular locations
#' @param nCol  number of columns of abundance levels
#' @param clustSize  average cluster size
#' @return assignProbsOut  assignments of each protein to compartments, and maximum (most common) assignment
protClusterStable <- function(geneProfileSummary, refLocProteins, refCol=1, nCol=9, clustSize=5, nReps=100, assignCutoff=0.5)   {

  # # # # # # # # #
  # check stability of kmeans by repeating it nReps times
  #  Since kmeans is a random process, this function will assess how
  #  stable the kmeans assignment process is
  #
  #  the file "refLocProteins" is a data frame with names "geneName" and "referenceCompartment"
  #  the variable "assignCutoff" gives the proportion (of nReps) of assignments needed to declare an assignment)
  # # # # # # # # #
  #nrowATdfm <- nrow(ATdfmTemp)
  nrowGeneProf <- nrow(geneProfileSummary)
  clMat <- data.frame(matrix("", ncol=nReps, nrow=nrowGeneProf))
  ntypes <- rep(NA, nrowGeneProf)
  purity.reps <- matrix(0, nrow=nrowGeneProf, ncol=8)   # 8 subcellular locations
  for (i in 1:nReps) {
    clProtsOrdAlltemp <- protCluster(geneProfileSummary=geneProfileSummary, refLocProteins=refLocProteins, refCol=refCol, nCol=nCol, clustSize=clustSize)
    ord.gene <- order(clProtsOrdAlltemp$clProtsOrd$geneName)    # alphabetical order by gene name
    clProtsOrdAll <- clProtsOrdAlltemp$clProtsOrd[ord.gene,]
    purityStatsOrd <- clProtsOrdAlltemp$purityStats[ord.gene,]
    purity.df.ord <- clProtsOrdAlltemp$purity.df[ord.gene,]

    cl.i <- purityStatsOrd$maxCat.vec
    clMat[,i] <- cl.i
    purity.reps <- purity.reps + as.matrix(purity.df.ord)
    print(i)
    if ((i %% 20) == 0) print(paste("replication ", i))
    }

  for (j in 1:nrowGeneProf) {
   #ntypes[j] <- length(unique(as.character(unlist(clMat[j,]))))
   table(as.character(unlist(clMat[1,])))
   }

  #table(ntypes)
  #table(as.character(unlist(clMat[3,])))

#  compartmentNames <- sort(as.character(unique(clProtsOrdAll$maxCat.vec)))
  #compartmentNames <- c("Cytosol",  "ER",   "Golgi", "Lysosome", "Mitochondrion",  "Nucleus",  "Peroxisome",  "PM", "none" )
  compartmentNames <- c(as.character(sort(unique(refLocProteins$referenceCompartment))), "unclassified")   # add "unclassified" to list of compartments
  compartmentCount <- function(xx, compartmentNames=compartmentNames) {
    # count numbers of proteins in each compartment for a particular set of nReps repetitions
    # "xx" is, for a particular protein, the list of assignments made by each of the nReps calls of the kmeans procedure
     xx <- as.character(unlist(xx))
     # xx <- clMat[3,]
     n.compartments <- length(compartmentNames)
     n.xx <- length(xx)
     n.c.j <- rep(0, n.compartments)
     for (i in 1:n.xx) {
       for (j in 1:n.compartments) {
         if (as.character(unlist(xx[i])) == compartmentNames[j]) n.c.j[j] <- n.c.j[j] + 1
         }
       }
     n.c.j
     }

  clMatTotals <- matrix(0, nrow=nrow(clMat), ncol=length(compartmentNames))
  maxAssign <- rep(NA, nrow(clMat))
  for (ii in 1:nrow(clMat)) {
    clMatTotals[ii,] <- compartmentCount(clMat[ii,], compartmentNames=compartmentNames)
    maxAssign[ii] <- maxLocAssign(clMatTotals[ii,], compartmentNames=compartmentNames, nReps=nReps, cutoff=assignCutoff)

    }
  clMatTotals.df <- data.frame(clMatTotals)
  names(clMatTotals.df)  <- compartmentNames
  clMatTotals2.df <- data.frame(clMatTotals.df, maxAssign)

  #geneNamesOrd <- order(geneProfileSummary$geneName)
  #geneProfileSummaryOrd <- geneProfileSummary[geneNamesOrd,]    # order by gene name, just to be sure

  OriginalClassif <- clProtsOrdAll$classif
  clMatTotals=data.frame(clMatTotals2.df, OriginalClassif)
  geneName <- geneProfileSummary$geneName

  row.names(geneProfileSummary) <- geneName
  row.names(clMatTotals) <- geneName
  row.names(purity.reps) <- geneName
  geneOrd <- order(geneName)    # fix, just in case

  geneProfileSummaryMclust <- list(geneProfileSummary=geneProfileSummary[geneOrd,],
                                         clMatTotals=clMatTotals[geneOrd,],
                                        purityMat=purity.reps[geneOrd,] )

  geneProfileSummaryMclust
  }

maxLocAssign <- function(matCounts, compartmentNames, nReps, cutoff=0.8) {
    # assign to category with propotion > cutoff; if none, call it "unclassified"
    # assignLocProps <- assignPropsAll[1,]
    if (cutoff <= 0.5) cutoff <- 0.5001
    exceedCutoff <- {matCounts/nReps > cutoff}
    if (sum(exceedCutoff) > 0) {
       indLoc <- which(exceedCutoff)
       protLoc <- compartmentNames[indLoc]
    }
    if (sum(exceedCutoff) == 0) protLoc <- "unclassified"
    protLoc
}

#' Compute precision and recall for assignments to subcellular compartments relative to reference
#'
#' @param assignVec list of assignments to subcellular compartments
#' @param compareVec list of reference ("true") assignments
#' @param nonName name given for genes with no assignment; must come alphabetically last
#'
#' @return cross-classification table of assigned vs comparison locations; then precision and recall

PrecisionRecallCompute <- function(assignVec, compareVec, noneName="unclassified") {
   # compare location assignments to a comparison (reference or "true") set of location assignments
   # unclassified genes must be labelled "unclassified", and must come alphabetically AFTER all assignment names
   #  Otherwise, there will be errors in the recall and precision calculations
  #cpaAssignVec.f <- assignVec
  #cpaTrueVec.f <- trueVec
  #assignTable <- table(cpaAssignVec.f,cpaTrueVec.f)
  assignVec <- as.character(assignVec)
  compareVec <- as.character(compareVec)
  compareVecOrig <- compareVec
  compareVec[compareVecOrig == noneName] <- NA
  assignVecOrig <- assignVec
  assignVec[assignVecOrig == noneName] <- NA
  assignTable <- table(assignVec, compareVec)
  rtotals <- margin.table(assignTable,1) # rows
  ctotals <- margin.table(assignTable,2) # cols
  assignTableR <- cbind(assignTable, rtotals)
  ctotalsAll <- c(ctotals, sum(ctotals))
  assignTableAll <- rbind(assignTableR, ctotalsAll)  # column names are "true" (reference) and row names "assigned"
  #write.csv(assignTableAll, file="assignTableTmtMs2Christo.csv")

  # now make tables including missing values in compareVec  (don't use for Recall and Precision)
  #  when compareVec is missing, it is NOT an error to classify it as having a certain location!
  #   (this is the point of the metamass routine)
  assignTableMiss <- table(assignVec, compareVec, useNA="ifany")
  rtotalsMiss <- margin.table(assignTableMiss,1) # rows
  ctotalsMiss <- margin.table(assignTableMiss,2) # cols
  assignTableRmiss <- cbind(assignTableMiss, rtotalsMiss)
  ctotalsAllmiss <- c(ctotalsMiss, sum(ctotalsMiss))
  assignTableAllmissing <- rbind(assignTableRmiss, ctotalsAllmiss)  # column names are "true" and row names "assigned"

  trueNames <- colnames(assignTable)
  compareNames <- rownames(assignTable)
  allLocNames <- union(trueNames, compareNames)
  #if (allLocNames[length(allLocNames)] != noneName) cat("error: classification names don't match\n")
  if (length(allLocNames[allLocNames != noneName]) != length(trueNames[trueNames != noneName])) {
    cat("error: true names and assign names don't match\n")
    }

  n.names <- length(trueNames)
  result <- NULL
  for (i in 1:n.names) {
    # i=2  # Cytosol
    TP <- as.numeric(assignTable[i,i])   # true positives
    Pos <- rtotals[i]        # all genes classified as positive
    FP <- Pos - TP           # false positives
    FN <-  ctotals[i] - TP   # false negatives (actually positvie, but classified as negative)
    Neg <- sum(ctotals) - Pos  # everything not positive
    TN <- Neg - FN           # true negatives
    Recall <- TP/(TP + FN)
    Precision <- TP/(TP + FP)
    result.i <- c(TP, FP, TN, FN, Pos, Neg, Recall, Precision)
    result <- cbind(result, result.i)
    }
  #result <- round(result, 1)
  colnames(result) <- trueNames
  rownames(result) <- c("TP", "FP", "TN", "FN", "Pos", "Neg", "Recall", "Precision")
  FalseNegCounts <- round(result[1:6,])
  RecallPrecision <- round(result[7:8,], 3)
  recallPrecisionAll <- rbind(FalseNegCounts, RecallPrecision)
  result <- list(assignTableAll=assignTableAll, recallPrecisionAll=recallPrecisionAll, assignTableAllmissing=assignTableAllmissing)
  result
   }

#result <- PrecisionRecallCompute(assignVec=cpaAssignVec.f, compareVec=cpaTrueVec.f)
