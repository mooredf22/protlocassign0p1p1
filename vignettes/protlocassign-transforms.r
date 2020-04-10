---
title: "Assignment of Proteins to Subcellular Locations Using protlocassign"
author: "Dirk Moore"
date: "`r Sys.Date()`"
#output: rmarkdown::html_vignette
output: word_document
vignette: >
  %\VignetteIndexEntry{protlocassign}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
Determining the locations of proteins in the cell is an important but complex problem. Differential centrifugation of cellular components into multiple fractions is a useful tool for doing this. This package takes data from centrifugation fractions for large numbers of proteins to determine their relative abundance among the fractions, and by extension to their relative abundance among subcelluar locations. The input data consisist of a file with one column, `geneName`, of gene/protein identifers, and multiple columns of the relative abundances of the proteins among the centrifugal fracions, with the abundance levels constrained to sum to 1. The fractions typically comprise the N, M, L1, L2, P, and S fractions, and may be supplemented by one to three additional "Nyc" fractons. Another key input to the program is a curated list of reference genes and their subcellular locations.This list is used to classify addtional genes to subcelluar locations.

This package implements a subcellular protein assignment procedure as "constrained proportional assignmet", or CPA (Jadot et al., 2016). This method first determines the profiles of abundance levels of genes in a curated reference set. Then, for each gene, it finds either a profile matching a single compartment or a linear combination of compartment profiles that match the profile of that gene. The relative weights of the linear combination in principle reflect the relative abundance of proteins for that gene among different compartments. The CPA method, unlike other previously described methods, can thus account for proteins that have multiple residences, and estimate the relative proportion among these residences.



## Tutorial on gene assignment

Let us illustrate with an example. Tannous et al. (2020) presented abundance levels of proteins among seven fractions: N, M, L1, L2, P, and S, and three additional fractions, "Nyc1", Nyc2", and "Nyc3", from Nycodenz fractions. Eight subcellular compartments were considered: nucleus, mitochondria, lysosomes, peroxisomes, golgi apparatus, plasma membrane, and cytosol. The CPA method assigns each of a large number of genes to one or more of these compartments, based on profiles from a set of reference genes. To run the program, the `prolocate` library must be installed. Also, the `BB` library is required; it may be downloaded from CRAN by typing `install.packages("BB") `.


```{r, echo=TRUE}
library(protlocassign)
dim(geneProfileSummaryTMTms2)
options(digits=3)
head(geneProfileSummaryTMTms2)
```
The data set gives relative fraction levels for 7893 genes. The last two columns give the number of spectra and the number of peptides (sequences) for each gene. (The spectra were averaged using a nested random effects model described in Jadot, 2016, to produce the given average for each gene.)

The reference genes, which are from Jadot et al (2016), may be accessed as follows:
```{r, echo=TRUE}
dim(refLocProteinsJadot)
refLocProteinsJadot[c(1,2,6,7,12,13,17,18,21,22,27,28,31,32,35,36),]
```
There are 39 genes in this reference set; two for each compartment are shown here.

To obtain profiles for the reference genes, use the function `cpaSetup`: 

```{r, echo=TRUE, eval=TRUE}
matLocRuse <- cpaSetup(geneProfileSummary=geneProfileSummaryTMTms2, refLocProteins=refLocProteinsJadot, n.channels=9)
matLocRuse

##referenceProfilePlot(refLocProteins = refLocProteinsJadot,
##  geneProfileSummary = geneProfileSummaryJadotExptA, matLocR = matLocR, n.channels=7)
```

To view them, use `referenceProfilePlot`:
```{r, echo=TRUE, eval=TRUE, fig.show='hold'}

#par(mfrow=c(4,3))
referenceProfilePlot(refLocProteins=refLocProteinsJadot, geneProfileSummary=geneProfileSummaryTMTms2,
                     matLocR=matLocRuse, n.channels=9)
```
Before finding the constrained proportional estimates, it is preferable to transform them to relative specific activities (RSA's), which are the amount of a given protein in a fraction divided by the amount of starting material of that protein (see Appendix section on adjusting protein levels). These adjustment factors are in the dataset "tmtMS2totProt", which is just a vector of length 9:

```{r, echo=TRUE, eval=TRUE, fig.show='hold'}
totProtUse <- tmtMS2totProt
totProtUse
```

The program "rsaDirect" adjusts the fraction levels measured from the machine to the amounts of protein in each fraction in the original sample. It is important to specify how many of the fractions are Nycodenz fractions (three in this case) so that the adjustment is just using the differential fracctions. We do this first with the reference profile levels in "matLocR":

```{r, echo=TRUE, eval=TRUE, fig.show='hold'}
matLocRrsa <- rsaDirect(geneProfileLevels=matLocRuse, nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
row.names(matLocRrsa) <- row.names(matLocRuse)
matLocRrsa
#round(matLocRrsa, digits=4)
```

We now do the same for the matrix of protein-level mean profiles:

```{r, echo=TRUE, eval=TRUE, fig.show='hold'}
geneProfileLevelsUse <- geneProfileSummaryTMTms2[,2:(1 + 9)]
geneProfileLevelsRSA <- rsaDirect(geneProfileLevels=geneProfileLevelsUse, 
                                   nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
dim(geneProfileLevelsRSA)
head(geneProfileLevelsRSA)
```

Now we put back the gene names in the first column, and Nspectra and Nseq in the last columns:

```{r, echo=TRUE, eval=TRUE, fig.show='hold'}
geneProfileSummaryRSA <- data.frame(geneProfileSummaryTMTms2[,1],geneProfileLevelsRSA, 
                                    geneProfileSummaryTMTms2[,11:12])
names(geneProfileSummaryRSA)[1] <- "geneName"
dim(geneProfileSummaryRSA)
head(geneProfileSummaryRSA)
```

Now we can plot the reference profiles on the RSA scale:

```{r, echo=TRUE, eval=TRUE, fig.show='hold'}
#par(mfrow=c(4,3))
referenceProfilePlot(refLocProteins=refLocProteinsJadot, geneProfileSummary=geneProfileSummaryRSA,
                     matLocR=matLocRrsa, n.channels=9)

```

Now we can run the CPA routine on the RSA-transformed levels; this may take several minutes to complete. The last column is a convergence indicator, with 1 indicating successful convergence. We remove this last column to get a matrix of gene names, and columns indicating the estimated proportional assignments of each gene among the eight subcellular locations.

```{r, echo=TRUE, eval=TRUE}
assignPropsOutTemp <- proLocAll(geneProfileSummary=geneProfileSummaryRSA, matLocR=matLocRrsa,
                            n.channels=9)
ncol.a <- ncol(assignPropsOutTemp)   # last column is a convergence indicator, which should be 1
assignPropsUse <- data.frame(assignPropsOutTemp[,-ncol.a] ) # drop this last column
head(assignPropsUse)
```



Next, load the original data, with spectra and peptides, and select the columns for gene name, peptide name, and the nine fractions ranging from N to Nyc.3:

```{r, echo=TRUE, eval=TRUE}
#finalListUse <- subset(tmtMS2orig, select={c(gene:peptide,N:Nyc.3)}, subset={gene %in% genesInclude})
finalListUse <- subset(tmtMS2orig, select={c(gene:peptide,N:Nyc.3)})
```

To convert to relative specific values, we convert this large dataframe in five smaller pieces, convert those pieces, and then reassemble:

```{r, echo=TRUE, eval=TRUE}
finalListRSAuse1 <- rsaDirect(geneProfileLevels=finalListUse[1:30000,3:11], 
                                   nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
finalListRSAuse2 <- rsaDirect(geneProfileLevels=finalListUse[30001:60000,3:11], 
                                   nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
finalListRSAuse3 <- rsaDirect(geneProfileLevels=finalListUse[60001:90000,3:11], 
                                   nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
finalListRSAuse4 <- rsaDirect(geneProfileLevels=finalListUse[90001:120000,3:11], 
                                   nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
finalListRSAuse5 <- rsaDirect(geneProfileLevels=finalListUse[120001:nrow(finalListUse),3:11], 
                                   nDiffFractions=6, nNycFractions=3, totProt=totProtUse)
finalListRSAuseT <- rbind(finalListRSAuse1, finalListRSAuse2, finalListRSAuse3, finalListRSAuse4, finalListRSAuse5)

finalListRSAuse <- data.frame(finalListUse[,1:2], finalListRSAuseT)

```
We can look at the profile of, for example, the gene "AADAC" by first finding the gene number in the dataset:
```{r, echo=TRUE, eval=TRUE}
protIndex("AADAC")
#protPlotfun(protPlot=93, assignProbsOut=assignProbsOut)
```
This function also accepts partial matching of the first few letters of a gene. For example, we can find the indices of the genes starting with "AAD""
```{r, echo=TRUE, eval=TRUE}
protIndex("AAD")
```
Now we plot the results for gene number AADAC, which has index number 93:

```{r, echo=TRUE, eval=TRUE, fig.width=7, fig.height=7}

protPlotfun(protPlot=93, geneProfileSummary=geneProfileSummaryRSA, finalList=finalListRSAuse, n.fractions=9,
                        Nspectra=T, matLocR=matLocRrsa, assignPropsMat=assignPropsUse)
```                        

The horizontal axis represents the nine fractions, which are N, M, L1, L2, P, S, Nyc.1, Nyc.2, and Nyc.3. In each of the eight plots, the red line is the average profile of the peptides, which are represented by the blue lines; thicker blue lines indicate larger numbers of spectra for a particular peptide. The dashed yellow-black lines show the expected profile for a protein entirely resident in the respective subcellular location. In this set of plots, we see that the CPA procedure assigns a 92\% residence proportion to endoplasmic reticulum. Visually, we see that the observed red profile closely matches the expected yellow-black profile for endoplasmic reticulum.

## Nearest genes in a data set
We may find the genes with profiles nearest to a given gene using the function "nearestGenes". Distance is computed as the Euclidean distance between profiles. To use the function, we first create a distance matrix for the genes in a list of mean profiles, such as "geneProfileSummaryTMTms2". In this dataset, the profiles are in columns 2 through 10.

```{r, echo=TRUE, eval=TRUE}
distUse <- dist(geneProfileSummaryTMTms2[,2:10], method="euclidean")
```
Then select the gene names:

```{r, echo=TRUE, eval=TRUE}
genesUse <- geneProfileSummaryTMTms2[,1]
```
Finally, provide a gene name. Here, for the gene "AADAC", we find the 10 nearest genes.
```{r, echo=TRUE, eval=TRUE}
protIndex("AAD")
nearestGenes("AADAC", n.nearest=10,  distGenes=distUse, geneNames=genesUse)
```

## Working with spectral level data
The spectral level data for tmtMS2 is in the dataset "tmtMS2orig". This dataset contains a list of 
spectra nested within peptides, which are in turn nested within proteins. Outliers may be identified as follows:

```
protsCombineCnew <- OutlierIdentify(protClass=tmtMS2orig, refCol=4, nCol=9, save.data=F)
```

The following function computes mean profiles using a random effects model with spectra nested within genes, 
after removing outliers:

```
gProfSummaryTMTms2 <- geneProfileEstimate(protsCombineCnew = protsCombineCnew, outlier="boxplot",
                                        refCol=4, nCol=9, eps=0.029885209, dataUse="tmtMS2", save.data=F)
head(gProfSummaryTMTms2)
```




## Christoforou example
Christoforou (2015) provided a set of profiles of 20 channels (from two separate runs) for over 5000 genes. Here are the constrained proportional assignments as computed by "proLocAll":


```{r, echo=TRUE, eval=TRUE}
geneProfileSummary=geneProfileSummaryChristoforou
refLocProteins=refLocProteinsChristoforou

```
Here are the the adjusted mean profiles. Here we only show the first four of the 20 columns.
```{r, echo=TRUE, eval=TRUE}
matLocRuse <- cpaSetup(geneProfileSummary=geneProfileSummaryChristoforou, refLocProteins=refLocProteinsChristoforou,
                       n.channels=20)
matLocRuse[,1:4]
```

Next we plot the reference profiles and compute the matrix of assigned proportions.

```{r, echo=TRUE, eval=TRUE}
referenceProfilePlot(refLocProteins=refLocProteinsChristoforou, geneProfileSummary=geneProfileSummaryChristoforou,
                     matLocR=matLocRuse, n.channels=20)

assignPropsOutTemp <- proLocAll(geneProfileSummary=geneProfileSummaryChristoforou, matLocR=matLocRuse,
                                n.channels=20)
ncol.a <- ncol(assignPropsOutTemp)   # last column is a convergence indicator, which should be 1
assignPropsUse <- data.frame(assignPropsOutTemp[,-ncol.a] ) # drop this last column
head(assignPropsUse)
```

# Appendix
## Outlier rejection and computation of mean profiles

The profiles of some compartments (e.g. nuclear and mitochondrial) show a high relative abundance level in the N or M fractions, respectively, and little in the others. Cytosolic proteins also are high in one fraction (S) and low in the others. Peroxisomes have high activity levels in the L1 and L2 fractions, and low levels in the others, including Nyc2. Lysosomes, like peroxisomes, have high level in the L1 and L2 fractions, but, unlike peroxisomes, high levels in Nyc2. The other compartments (ER, Golgi, and PM) have more complex profiles. Once we have obtained reference profiles for each of the eight compartments, we score each new protein as the optimal combination of the reference proteins, thereby obtaining estimates of the proportionate contributions of each reference profile to that of the protein being considered.


The first step is an outlier screen, for each protein, of the abundance levels for each fraction of all of the spectra. Abundance levels $p$ are first log2 transformed via $y = log_2(p + \delta)$, where $\delta$ is an estimate of the background "noise" in the abundance estimates. Here we used $\delta = 10^{(-5)}$. Then boxplot outliers are identified as values more than three times the interquartile range beyond the first or third quartile. We discard any spectrum for which any observation is thereby classified as an outlier. This process is repeated for all proteins.


## Constrained proportional assignment
Once outlier spectra have been removed, the next step is to determine a mean profile for each protein. If a protein has at least 4 spectra and at least 3 different sequences (peptides), we ordinarily have enough data to fit a random effects model. This model computes, for each EF in each fraction of this protein, a weighted average and standard error of the measures, accounting for the fact that spectra are nested within sequences. (This computation is carried out using the “lmer” function in the “lme4” R package.) . The result is essentially the mean and standard error of the observations, with an adjustment for the nested structure of the data. This procedure prevents a sequence with a very large number of spectra from dominating the estimates of the mean and standard error.  Occasionally the “lmer” program will fail to converge for a particular abundance. If that happens or if the condition of having at least 4 spectra and 3 different sequences is not met, we compute the simple mean and standard error of the (log2 transformed) for that abundance level. If there are fewer than 3 spectra, we cannot compute a standard error, and thus only report the mean.  When this process has been completed, every protein has a mean profile $\underset{\sim}x = (x_1, x_2, ..., x_q)$  of abundance levels in the $q$ fractions N, M, L1, L2, P, S, and Nyc2.  For proteins with at least 2 peptides and three spectra, each component of the profile has an estimate of its standard deviation. 

Next, the profiles of the reference proteins are selected. For each of the eight compartments (Cytosol, ER, Golgi, Lysosome, Mitochondria, Nucleus, Peroxisome, and PM), the reference protein relative abundance levels are averaged to form eight compartment profiles $\underset{\sim}s_1, \underset{\sim}s_2, ..., \underset{\sim}s_8$, where each vector $\underset{\sim}s_j$ is a $q \times 1$ vector of mean levels of the $q=7$ or $q=9$ channels. 

Finally, for each protein, we find estimated proportions $\hat{p}_1, \hat{p}_2, ..., \hat{p}_8$ so that
$$\underset{\sim}y = \hat{p}_1 \underset{\sim}s_1 + ...+ \hat{p}_8 \underset{\sim}s_8$$
is as close as possible to the observed value $\underset{\sim}x$, subject to the constraints

$$ 0 \leq p_j \leq 1$$

for all $j$, and



$$ \sum_{j=1}^{8} \hat{p}_j = 1 $$
where "close" is defined by minimizing the sum of squares of the differences

$$ Q=\sum_{i=1}^{q}(y_i - x_i)^2 $$
Thus, we may view the proportions $\hat{p}_j$ as proportional allocations of the eight standard profiles to form $\underset{\sim}y$, which is as close as possible to the observed $\underset{\sim}x$ for this particular protein. This constrained optimization is carried out using the “spg” function in the R package “BB” to compute assignment probabilities for each profile for each organelle. 

## Adjustment of protein abundance levels

The input data, geneProfileSummary, consists of a list of gene names (first column) and the relative abundance of the corresponding proteins in the next columns. Typically, there will be six differential fraction columns (N, M, L1, L2, P, and S) and one to three additional Nycodenz fractions (which have been extracted from the L1 differential fraction). These quantities represent the relative amounts of a protein in each fraction. The amounts of protein extracted from each fraction, totProtein, differ dramatically from one fraction to another. Notably, the L1 and L2 fractions, which are heavily enriched in lysosomal and peroxixomal proteins, contain much smaller amounts of protein than the other differential fractions. Also, the samples used for the Nycodenz fractionization are taken from the L1 fraction, so they also represent very small proportions of the total sample. For technical reasons, when samples of material are prepared for analysis, equal amounts are selected for each channel, regardless of the amount of protein in the fractions. As a result, when the data are analyzed, adjustments need to be made in order to obtain meaningful comparisons.

The first transformation involves finding the amount of a given protein in a fraction divided by the amount of protein in the starting material. To see how this works, let us consider the Jadot reference protein profiles

```{r, echo=TRUE, eval=TRUE}
matLocR <- cpaSetup(geneProfileSummary=geneProfileSummaryTMTms2, refLocProteins=refLocProteinsJadot, n.channels=9)
round(matLocR, digits=4)

##referenceProfilePlot(refLocProteins = refLocProteinsJadot,
##  geneProfileSummary = geneProfileSummaryJadotExptA, matLocR = matLocR, n.channels=7)
```

Here, each row represents the profile for a protein resident solely in a particular compartment. The amount of starting material, tmtMS2totProt, in each fraction is as follows:

```{r, echo=TRUE, eval=TRUE}
#totProt=c(46.044776, 48.955954, 1.384083, 1.566324, 24.045584, 58.181818, 0.0368564, 0.0684596, 1.27301587)
totProtUse <- tmtMS2totProt
totProtUse

totProtdf <- data.frame(t(matrix(totProtUse)))
names(totProtdf) <- colnames(matLocR)
totProtdf
```
The total amount of starting material is the sum of the amounts given in the first six (differential) fractions (N, M, L1, L2, P, and S).
```{r, echo=TRUE, eval=TRUE}
sum(totProtUse[1:6])
```
(The last three "Nyc" columns were derived from L1.)

To compute the amount of protein in each fraction, we multiply each column by the amount of starting material in each fraction. The function abundanceTransform does this:

```{r, echo=TRUE, eval=TRUE}
#totProt=c(46.044776, 48.955954, 1.384083, 1.566324, 24.045584, 58.181818, 0.0368564, 0.0684596, 1.27301587)

protAbund <- abundanceTransform(matLocR,6,3, totProt=totProtUse)
round(protAbund$amtProtFrac, digits=4)
round(protAbund$relAmtProtFrac, digits=4)
```

The first component of the result, amtProtFrac, is the amount of protein in each fraction for each row. The second component, relAmtProtFrac,
is the amount of given protein in fraction / amount of given protein in starting material. For example, for a cytosolic protein, the amount of protein in each fraction is given by

```{r, echo=TRUE, eval=TRUE}
matLocR[,1]*totProtUse[1]
protAbund$amtProtFrac[1,]
```
This is the first column of amtProtFrac.

The first row of relAmtProtFrac, to continue with the example, is the first row of amtProtFrac (amount of cystolic protein in each fraction) divided by the total amount of cystolic protein, which is the sum of the first row of amtProtFrac:

```{r, echo=TRUE, eval=TRUE}
protAbund$amtProtFrac[1,] / sum(protAbund$amtProtFrac[1,])
```

The relative amount of a cytosolic protein is given by the amounts of protein in the first row (a cytosolic protein) divided by the amounts in all six differential fractions (the first six elements of the first row):
```{r, echo=TRUE, eval=TRUE}

protAbund$amtProtFrac[1,]/sum(protAbund$amtProtFrac[1,])
protAbund$relAmtProtFrac[1,]
```
This is the first row of relAmtProtFrac.

The values in relAmtProtFrac represent the amount of protein per fraction, normalized to same amount of protein in the starting material (before mixing). In theory, these values correspond approximately to the relative amounts of protein in compartments within a cell.

Finally, we compute the relative specific activity (RSA):

To get this, we first find Difp, the total protein in the six differential fractions (nDiffFractions = 6), and then the proportions propFrac of protein in the differential fractions:
```{r, echo=TRUE, eval=TRUE}
rsaResult <- RSAtransform(protAbund$relAmtProtFrac, totProt=totProtUse)
round(rsaResult$rsa, digits=4)
```

```{r, echo=TRUE, eval=TRUE}
Difp <- sum(totProtUse[1:6])   # total protein in the differential fractions
Difp
propFrac <- totProtUse/Difp  # proportion of protein in the differential fractions
propFrac
```

For example, for the N fraction (first column), we have, which matches the first column of rsa:
```{r, echo=TRUE, eval=TRUE}
protAbund$relAmtProtFrac[,1]/propFrac[1]
rsaResult$rsa[,1]
```

Finally, if we normalize the rsa matrix so that rows sum to one, we get original input data, matLocR
```{r, echo=TRUE, eval=TRUE}
round(rsaResult$rsaFractions, digits=4)
round(matLocR, digits=4)
```

If we just need to obtain the RSA transformed data directly, just do this:
```{r, echo=TRUE, eval=TRUE}
matLocRrsa <- rsaDirect(geneProfileLevels=matLocR, nDiffFractions=6, nNycFractions=3, 
        totProt=totProtUse)

round(matLocRrsa, digits=4)
```

## Simulating proteins resident in multiple subcellular locations

We may simulate data from proteins with multiple residences using the "proteinMix" function. For example, to simulate date from proteins resident in a range of proportions in cytosol and lysosomes, we do this:
```{r, echo=TRUE, eval=TRUE}
#totProt=c(46.044776, 48.955954, 1.384083, 1.566324, 24.045584, 58.181818, 0.0368564, 0.0684596, 1.27301587)

protAbund <- abundanceTransform(matLocR, nDiffFractions=6, nNycFractions=3 , totProt=totProtUse)
relAmtProtFrac <- protAbund$relAmtProtFrac
mixCytoLyso <- proteinMix(relAmtProtFrac, 1, 4)
round(mixCytoLyso, digits=4)
```

Then we can test the CPA algorithm by first converting this mixture data to RSA's:
```{r, echo=TRUE, eval=TRUE}

mixCytoLysoRSA <- RSAtransform(relAmtProtFrac=mixCytoLyso, nDiffFractions=6, nNycFractions=3, totProt=totProtUse)$rsa
  
mixCytoLysoRelAmt <- data.frame(rownames(mixCytoLysoRSA), mixCytoLysoRSA)
names(mixCytoLysoRelAmt)[1] <- "geneName"
round(mixCytoLysoRelAmt[,-1], digits=4)

#matLocRrsa <- rsaDirect(geneProfileLevels=matLocR, totProt=totProt)
matLocRrsa <- rsaDirect(geneProfileLevels=matLocR, nDiffFractions=6, nNycFractions=3, totProt=totProtUse, maxRSA=25)
rownames(matLocRrsa) <- rownames(matLocR)
```
Finally, we fit the CPA algorithm to this RSA-transformed, simulated data:

```{r, echo=TRUE, eval=TRUE}
mixCytoLysoPropT <- proLocAll(geneProfileSummary=mixCytoLysoRelAmt, matLocR=matLocRrsa,
                            n.channels=9)
mixCytoLysoProp <- mixCytoLysoPropT[,-c(1,10)]
row.names(mixCytoLysoProp) <- mixCytoLysoPropT[,1]
round(mixCytoLysoProp, digits=4)

```
The estimated proportions correspond closely to the the proportions used in the simulation.

## References

Christoforou A, Mulvey, CM, Breckels LM, Geladaki A, Hurrell T, Hayward PC, Naake T, Gatto L, Viner R, Arias AM, and Lilley KS (2016), "A draft map of the mouse pluripotent stem cell spatial proteome" Nature Communications 7, 8992. DOI: 10.1038/ncomms9992

Jadot M, Boonen M, Thirion J, Wang N, Xing J, Zhao C, Tannous A, Qian M, Zheng H, Everett JK, Moore DF, Sleat DE, Lobel P (2016) Accounting for protein subcellular localization: a compartmental map of the rat liver proteome. Molecular and Cellular Proteomics 16, 194-212. doi:10.1074/mcp.M116.064527  PMCID: PMC5294208

Tannous A, Boonen M, Zheng H, Zhao C, Germain C, Moore D, Sleat D, Jadot M, Lobel P. Comparative Analysis of Quantitative Mass Spectrometric Methods for Subcellular Proteomics. Journal of Proteome Research 2020, in press

