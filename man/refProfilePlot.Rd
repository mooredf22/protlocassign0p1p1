% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refProfilePlot.r
\name{refProfilePlot}
\alias{refProfilePlot}
\title{Plot profiles of reference proteins}
\usage{
refProfilePlot(
  refLoc,
  markerList = markerListJadot,
  protProfileSummary = protProfileSummaryTMTms2,
  markerProfiles = markerProfiles,
  refProtPlot = NULL
)
}
\arguments{
\item{refLoc}{the name of the reference subcellular compartment to plot}

\item{markerList}{List of reference proteins}

\item{protProfileSummary}{data frame of protein names (as row names) and relative abundance levels.}

\item{markerProfiles}{A matrix markerProfiles giving the abundance level profiles of the subcellular locations}
}
\description{
This function plots profiles of reference proteins and also the average profile for each compartment
}
