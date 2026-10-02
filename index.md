# lblvm

*Please feel free to create a* [GitHub
issue](https://github.com/SeewooLi/lblvm/issues) *for bug reports or
potential improvements.*

## About the Framework

The **logit-beta latent variable modeling (lblvm) framework** can be
viewed alongside two well-established latent variable modeling
traditions: **factor analysis (FA)** and **item response theory (IRT)**.
Conventional FA models are typically formulated with an identity link
and the normally distribution, whereas IRT generally uses a logit link
(or a probit link) with the Bernoulli or the categorical distribution.

The underlying premise of the logit-beta framework is parallel to that
of FA and IRT: observed responses can be interpreted as manifestations
of one or more **unobserved continuous latent variables**. The framework
is conceptually similar to FA in modeling responses with a continuous
probability distribution, while sharing with IRT the property of
**mean-dependent conditional variance**. Unlike the homoscedastic normal
distribution typically used in FA, the variance of a beta distribution
depends on its mean. As a result, the intercept plays a more prominent
role, similar to its role in IRT.

The framework accommodates both **bounded-continuous and ordered
categorical responses**. Bounded-continuous responses can be naturally
modeled using the beta distribution ([Li & Shin,
2025](https://doi.org/10.1017/psy.2025.10044)). For ordered categorical
responses, the **discretized response model (DRM)** assumes that
observed categories arise through the discretization of an underlying
bounded-continuous response process, with category probabilities
determined by intervals on the bounded continuum (Li & Jeon, in review).

## Package Overview

`lblvm` provides tools for the logit-beta framework, including parameter
estimation, model comparison, and visualization.

## Installation

You can install the development version of `lblvm` from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("SeewooLi/lblvm")
```
