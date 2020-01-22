\name{tmtMS2totProt}
\alias{tmtMS2totProt}
\docType{data}
\title{tmtMS2totProt
%%   ~~ data name/kind ... ~~
}
\description{
Adjustment factors for the tmtMS2 data, to adjust the relative protein amounts as measured
by the mass spectrometer to the amounts of protein in each channel in the original sample.
This is a vector with nine values: there is one value for each of the six fraction channels and the three Nycodenz channels.
Note that the sum of the three Nycodenz channels (Nyc.1, Nyc.2, Nyc.3) equals the amount in L1.
These values are used to adjust the levels in the data set "geneProfileSummaryTMTms2" to re-create the levels in the original samples.
%%  ~~ A concise (1-5 lines) description of the dataset.
}
\usage{data("tmtMS2totProt")}
\format{
  The format is:
 names: "N" "M" "L1" "L2" "P" "S" "Nyc.1" "Nyc.2" "Nyc.3""
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
A Tannous et al. (2020) submitted
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
tmtMS2totProt
}
\keyword{datasets}
