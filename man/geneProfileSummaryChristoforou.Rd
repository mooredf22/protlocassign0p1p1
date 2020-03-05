\name{geneProfileSummaryChristoforou}
\alias{geneProfileSummaryChristoforou}
\docType{data}
\title{
geneProfileSummaryChristoforouExpt1
%%   ~~ data name/kind ... ~~
}
\description{This is a list of 5032 genes and their relative expression levels in Experiments 1 and 2, described in Christoforou et al (2016). The first column is the gene name, and the next two sets of ten columns are the channel values for experiments 1 and 2, respectively

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("geneProfileSummaryChristoforouExpt1")}
\format{
  A data frame with 5032 observations on the following 12 variables.
  \describe{
    \item{\code{geneName}}{a character vector representing the gene name}
    \item{\code{X126}}{a numeric vector}
    \item{\code{X127N}}{a numeric vector}
    \item{\code{X127C}}{a numeric vector}
    \item{\code{X128N}}{a numeric vector}
    \item{\code{X128C}}{a numeric vector}
    \item{\code{X129N}}{a numeric vector}
    \item{\code{X129C}}{a numeric vector}
    \item{\code{X130N}}{a numeric vector}
    \item{\code{X130C}}{a numeric vector}
    \item{\code{X126.1}}{a numeric vector}
    \item{\code{X127N.1}}{a numeric vector}
    \item{\code{X127C.1}}{a numeric vector}
    \item{\code{X128N.1}}{a numeric vector}
    \item{\code{X128C.1}}{a numeric vector}
    \item{\code{X129N.1}}{a numeric vector}
    \item{\code{X129C.1}}{a numeric vector}
    \item{\code{X130N.1}}{a numeric vector}
    \item{\code{X130C.1}}{a numeric vector}
    \item{\code{Nspectra}}{a numeric vector indicating number of spectra}
    \item{\code{Nseq}}{a numeric vector indicating number of peptides (sequences)}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Chrisoforou A, Mulvey CM, Breckels LM, Geladaki A, Hurrell T, Hayward PC, Naake T, Gatto L, Viner R, Arias AM, Lilley KS (2016) A draft map of the mouse pluripotent stem cell spatial proteome. Nature Communications
 https://doi.org/10.1038/ncomms9992

%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

head(geneProfileSummaryChristoforou)
}
\keyword{datasets}
