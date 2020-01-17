\name{tmtMS2orig}
\alias{tmtMS2orig}
\docType{data}
\title{
tmtMS2orig
%%   ~~ data name/kind ... ~~
}
\description{
This is the full dataset for TMT-MS2 (A Tannous et al., 2020), including genes, peptides, and spectra
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("tmtMS2orig")}
\format{
  A data frame with 131699 observations on the following 13 variables.
  \describe{
    \item{\code{gene}}{a character vector}
    \item{\code{peptide}}{a character vector}
    \item{\code{N}}{a numeric vector}
    \item{\code{M}}{a numeric vector}
    \item{\code{L1}}{a numeric vector}
    \item{\code{L2}}{a numeric vector}
    \item{\code{P}}{a numeric vector}
    \item{\code{S}}{a numeric vector}
    \item{\code{Nyc.1}}{a numeric vector}
    \item{\code{Nyc.2}}{a numeric vector}
    \item{\code{Nyc.3}}{a numeric vector}
    \item{\code{geneId}}{a numeric vector}
    \item{\code{geneSeqId}}{a numeric vector}
  }
}
\details{
geneID is a unique numeric identifier for gene, and geneSeqId is a unique numeric identifier for peptide
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
A Tannous et al. (2020) submitted manuscript
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(tmtMS2orig)
str(tmtMS2orig)
}
\keyword{datasets}
