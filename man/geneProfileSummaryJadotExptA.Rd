\name{geneProfileSummaryJadotExptA}
\alias{geneProfileSummaryJadotExptA}
\docType{data}
\title{geneProfileSummaryJadotExptA
%%   ~~ data name/kind ... ~~
}
\description{
A data frame with 8071 observations on relative protein levels in 7 channels from Jadot Experiment A.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("geneProfileSummaryJadotExptA")}
\format{
  A data frame with 8071 observations on the following 10 variables.
  \describe{
    \item{\code{geneName}}{a character vector}
    \item{\code{N}}{a numeric vector}
    \item{\code{M}}{a numeric vector}
    \item{\code{L1}}{a numeric vector}
    \item{\code{L2}}{a numeric vector}
    \item{\code{P}}{a numeric vector}
    \item{\code{S}}{a numeric vector}
    \item{\code{Nyc2}}{a numeric vector}
    \item{\code{Nspectra}}{a numeric vector}
    \item{\code{Nseq}}{a numeric vector}
  }
}
\details{Nspectra is a count of the total number of spectra for each gene, and Nseq is the number of sequences (peptides) for that gene.

%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Experimental details may be found in Jadot M, Boonen M, Thirion J, Wang N, Xing J, Zhao C, Tannous A, Qian M, Zheng H, Everett JK, Moore DF, Sleat DE, and Lobel P (2016) Accounting for protein subcellular localization: a compartmental map of the rat liver proteome. Molecular and Cellular Proteomics 16, 194-212. doi:10.1074/mcp.M116.064527 PMCID: PMC5294208
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(geneProfileSummaryJadotExptA)
## maybe str(geneProfileSummaryJadotExptA) ; plot(geneProfileSummaryJadotExptA) ...
}
\keyword{datasets}
