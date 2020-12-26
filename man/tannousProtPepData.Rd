\name{tannousProtPepData}
\alias{tannousProtPepData}
\docType{data}
\title{tannousProtPepData
%%   ~~ data name/kind ... ~~
}
\description{
This data set contains profiles of all proteins and peptides. The format is first a
protein profile, followed by all derivative peptide profiles
%%  ~~ A concise (1-5 lines) description of the dataset.
}
\usage{data("tannousSpectral")}
\format{
  The format is:
 names: "prot" "peptide" "N" "M" "L1" "L2" "P" "S" "Nyc.1" "Nyc.2" "Nyc.3""
}
\details{
Row names are either protein or peptide. For the protein record, the protein name is
repeated in the peptide column. The outlier column identifies outlier peptides.
Other columns relevant to the protein are also included, as in the data set tannousSpectra.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
A Tannous et al. (2020)
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
str(tannousProtPepData)
}
\keyword{datasets}
