# Bayes Factor from Laplace Approximation

Bayes Factor from Laplace Approximation

## Usage

``` r
BF(fit0, fit1)
```

## Arguments

- fit0:

  Model 0

- fit1:

  Model 1

## Value

A list containing the following objects:

- log_marginal_M0:

  Logarithm of the marginal probability of Model 0.

- log_marginal_M1:

  Logarithm of the marginal probability of Model 1.

- logBF:

  Logarithm of the Bayes factor

- BF:

  The Bayes factor

## Details

Let \\D\\ denote the observed data, \\\tau\\ the parameter vector, and
\\M_k\\ a model or hypothesis indexed by \\k\\. The marginal probability
of the data under \\M_k\\ is \$\$ I_k = P(D \mid M_k) = \int P(D \mid
\tau, M_k)\pi(\tau \mid M_k)\\d\tau, \$\$ where \\\pi(\tau \mid M_k)\\
is the prior density of \\\tau\\. The marginal probability is
approximated using a Laplace approximation: \$\$ I_k \approx
(2\pi)^{d/2} \lvert\hat{\Sigma}\rvert^{1/2} P(D \mid \hat{\tau}, M_k)
\pi(\hat{\tau} \mid M_k), \$\$ where \\d\\ is the dimension of the
parameter vector, \\\hat{\tau}\\ is the mode of the posterior
distribution, and \\\hat{\Sigma}^{-1}\\ is the Fisher information
matrix. For example, `log_marginal_M0` denotes \\\log I_0\\, the log
marginal probability of the data under Model 0.
