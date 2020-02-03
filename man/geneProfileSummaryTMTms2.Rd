\name{geneProfileSummaryTMTms2}
\alias{geneProfileSummaryTMTms2}
\docType{data}
\title{Protein mean relative abundance levels from the TMT-MS2 data
%%   ~~ data name/kind ... ~~
}
\description{Each of the 7893 genes (proteins) has relative abundance levels in the channels N, M, L1, L2, P, S, Nyc1, Nyc2, and Nyc3, constrained to sum to 1. Details may be found in A. Tannous et al. (2020). These levels refer to the relative levels measured by the mass spectrometer, averaged across spectra and peptides for each gene. They are not adjusted for protein amounts derived from the original samples.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("geneProfileSummaryTMTms2")}
\format{
  A data frame with 7893 observations on the following 12 variables.
  \describe{
    \item{\code{geneName}}{a factor with levels for each of the 7893 genes }
    \item{\code{N}}{a numeric vector}
    \item{\code{M}}{a numeric vector}
    \item{\code{L1}}{a numeric vector}
    \item{\code{L2}}{a numeric vector}
    \item{\code{P}}{a numeric vector}
    \item{\code{S}}{a numeric vector}
    \item{\code{Nyc.1}}{a numeric vector}
    \item{\code{Nyc.2}}{a numeric vector}
    \item{\code{Nyc.3}}{a numeric vector}
    \item{\code{Nspectra}}{a numeric vector}
    \item{\code{Nseq}}{a numeric vector}
    }
  }

\details{
Nspectra is a count of the total number of spectra for each gene, and Nseq is the number of sequences (peptides) for that gene.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Tannous A, Boonen M, Zheng H, Zhao C, Germain C, Moore D, Sleat D, Jadot M, Lobel P. Comparative Analysis of Quantitative Mass Spectrometric Methods for Subcellular Proteomics. Journal of Proteome Research 2020, in press
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Experimental details may be found in Jadot M, Boonen M, Thirion J, Wang N, Xing J, Zhao C, Tannous A, Qian M, Zheng H, Everett JK, Moore DF, Sleat DE, and Lobel P (2016) Accounting for protein subcellular localization: a compartmental map of the rat liver proteome. Molecular and Cellular Proteomics 16, 194-212.
doi:10.1074/mcp.M116.064527  PMCID: PMC5294208
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(geneProfileSummaryTMTms2)
## maybe str(geneProfileSummaryTMTms2) ; plot(geneProfileSummaryTMTms2) ...
}
\keyword{datasets}
