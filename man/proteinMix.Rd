% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSAtransform.r
\name{proteinMix}
\alias{proteinMix}
\title{compute a mixture of proteins in two compartments}
\usage{
proteinMix(AcupRef, Loc1, Loc2, increment = 0.1)
}
\arguments{
\item{AcupRef}{amount of given protein in fraction / amount of given protein in starting material}

\item{Loc1}{row number of one compartment}

\item{Loc2}{row number of other compartment}

\item{increment}{fraction increment from 0 to 1}
}
\value{
mixAmount relative amounts of proteins in the fractions
}
\description{
compute a mixture of proteins in two compartments
}
