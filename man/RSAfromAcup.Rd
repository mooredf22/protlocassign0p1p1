% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSAtransform.r
\name{RSAfromAcup}
\alias{RSAfromAcup}
\title{Compute relative specific amount from relative amounts in protein fractions}
\usage{
RSAfromAcup(Acup, NstartMaterialFractions = 6, totProt = NULL)
}
\arguments{
\item{Acup}{amount of given protein in fraction / amount of that given protein in starting material}

\item{NstartMaterialFractions}{Number of differential fractions, typically 6, for N, M, L1, L2, P, and S}

\item{totProt}{Total protein counts in each of the differential and nycodenz fractions; this is necessary to compute RSA's}
}
\value{
rsa: relative specific amount
}
\description{
Compute relative specific amount from relative amounts in protein fractions
}
