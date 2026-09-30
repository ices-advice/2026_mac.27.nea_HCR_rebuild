Shortcut MSE for NEA mackerel (mac.27.nea) 2026
================

Summary: see
[`mac.27.nea_OM_2026.md`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/mac.27.nea_OM_2026.md)
for the operating model and
[`mac.27.nea_MP_2026.md`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/mac.27.nea_MP_2026.md)
for the results.

## Introduction

This repository contains the data and code for a **shortcut** MSE for:

- Mackerel (*Scomber scombrus*) in subareas 1-8 and 14 and in divisions
  9.a, 12.a, and 12.b (Northeast Atlantic and adjacent waters) -
  mac.27.nea

This shortcut MSE was developed in response to a request to ICES “To
provide a technical service regarding the impacts of different
rebuilding Harvest Control Rules (HCRs) for Altantic mackerel
(mac.27.nea)”.

*This work is only intended as an interim analysis until a full MSE is
completed.*

The approach of generating the operating model (OM) is the standard
approach used in ICES MSEs, i.e. sampling from the variance-covariance
matrix from the SAM model fit. This is the same approach as used in
recent ICES MSEs, e.g. WKNSMSE 2019
(<https://doi.org/10.17895/ices.pub.5090>), WKBPLAICE 2025
(<https://doi.org/10.17895/ices.pub.28400255>), WKMSEHERRING 2025
(<https://doi.org/10.17895/ices.pub.28868798>). Details of the approach
are e.g. described in this WKBPLAICE 2025 working document:
<https://github.com/shfischer/WKBPLAICE2024_ple.27.7e_MSE/blob/WKBPLAICE2024/WKBPLAICE2024_ple.27.7e_OM.pdf>

The analysis is based on FLR’s mse package.

## Repository structure

- `funs_*`: Function libraries, defining the functions used in the other
  scripts

  - [`funs_MSE.R`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/funs_MSE.R):
    functions (MSE modules) used in the MSE projection

  - [`funs_OM.R`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/funs_OM.R):
    functions used for conditioning the OMs

- `mac.27.nea_*.md`: R markdown documents with the analysis (code,
  plots, descriptions). The files are provided as `.Rmd` (R code) and
  `.md` (rendered documents for showing on GitHub)

  - [`mac.27.nea_OM_2026.md`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/mac.27.nea_OM_2026.md):
    document for generating the OM (update, based on WGWIDE 2026
    assessment)

  - [`mac.27.nea_MP_2026.md`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/mac.27.nea_MP_2026.md):
    document for running the MSE and analysing the results (update,
    based on WGWIDE 2026 assessment)

  - [`mac.27.nea_OM.md`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/mac.27.nea_OM.md):
    document for generating the OM (based on WGWIDE 2025 assessment)

  - [`mac.27.nea_MP.md`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/mac.27.nea_MP.md):
    document for running the MSE and analysing the results (based on
    WGWIDE 2025 assessment)

- `renv.lock`: renv file specifying the version of R and R packages (see
  below)

[`input/`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/tree/main/input):
This directory contains all files required for generating the OM

- [`input/fit.rds`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/input/fit.rds):
  SAM model fit from the WGWIDE 2025 stock assessment

- [`input/fit_2026.rds`](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/input/fit_2026.rds):
  SAM model fit from the WGWIDE 2026 stock assessment

## R, R packages and version info

The MSE simulations were run with R:

``` r
> sessionInfo()
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)
```

The package versions and their dependencies are recorded with the R
package [renv](https://rstudio.github.io/renv/) and stored in the file
[renv.lock](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/renv.lock).
The exact package version can be restored by cloning this repository,
navigating into this folder in R (or setting up a project), installing
the renv package

``` r
install.packages("renv")
```

and calling

``` r
renv::restore()
```

See [renv](https://rstudio.github.io/renv/) and the package
documentation for details.

The framework is based on the Fisheries Library in R (FLR) framework and
uses the [FLR packages](https://flr-project.org/)
[`FLCore`,](https://github.com/flr/FLCore)
[`FLasher`](https://github.com/flr/FLasher),
[`FLBRP`](https://github.com/flr/FLBRP),
[`ggplotFL`](https://github.com/flr/ggplotFL),
[`mse`](https://github.com/flr/mse), and
[`FLfse`](https://github.com/shfischer/FLfse). See
[renv.lock](https://github.com/ices-taf/2026_mac.27.nea_HCR_rebuild/blob/main/renv.lock)
for version details and sources.

Also, the R package
[`stockassessment`](https://github.com/fishfollower/SAM)is used.
