# protlocassign: A method for proportionally assigning protein residence among subcellular organelles

## Introduction
This package implements constrained proportional assignment to assign proteins proportionately to their subcellular organelle residences. It uses high througput mass spectrometry data from serial centrifugal fractionation of tissue samples to determine the relative distribution of proteins in the six deDuve fractions N, M, L1, L2, P, and S, and in three Nycodenz fractions (derived from L1), namely Nyc1, Nyc2, and Nyc3. The package identifies outlier spectra and computes weighted means across peptides using random effects models. Then, using reference profiles, it estimates the distribution of each protein across the main compartments by obtaining optimal mixtures of the reference profiles. 

## Installation
Start by installing the devtools package from CRAN, by typing

```
install.packages("devtools")
```

Then install the protlocassign package by typing

```
devtools::install_github("protlocassign0p1p1")
library(protlocassign0p1p1)
```

This will make the programs and datasets available.Instructions for use of the package are in the vignette.
