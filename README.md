
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RiskMap

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/RiskMap)](https://cran.r-project.org/package=RiskMap)
[![R-CMD-check](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/claudiofronterre/RiskMap/actions/workflows/R-CMD-check.yaml)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/RiskMap)](https://cran.r-project.org/package=RiskMap)
[![Total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/RiskMap)](https://cran.r-project.org/package=RiskMap)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

## Overview

`RiskMap` provides tools for model-based geostatistical analysis of
continuous, binomial and Poisson outcomes.

- Fit spatial and spatio-temporal Gaussian process models.
- Generate predictive surfaces and target summaries.
- Run simulation-based diagnostics and validation workflows.

The methodology is described in *Model-based Geostatistics for Global
Public Health* by Diggle and Giorgi.

## Start Here: MBG-R Book

For a full applied guide to using `RiskMap` in real public health
workflows, see the online book by **Emanuele Giorgi and Claudio
Fronterre**:

**Model-based geostatistics for global public health using R**
<https://www.mbgr.org/>

## Installation

Install the stable version from CRAN:

``` r
install.packages("RiskMap")
```

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("claudiofronterre/RiskMap")
```

## Quickstart

A minimal linear Gaussian geostatistical model of the form
$Y(x) = \beta_0 + S(x)$, where $S(x)$ is a spatial Gaussian process can
be fitted with:

``` r
library(RiskMap)

data(italy_sim)

fit <- glgpm(
  formula = y ~ gp(x1, x2),
  data = italy_sim,
  family = "gaussian",
  crs = 32634,
  messages = FALSE
)

summary(fit)
```

## Learn More

- MBG-R book: <https://www.mbgr.org/>
- Package website: <https://claudiofronterre.github.io/RiskMap/>
- Issue tracker: <https://github.com/claudiofronterre/RiskMap/issues>
