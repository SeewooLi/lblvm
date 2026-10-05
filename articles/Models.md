# Models

This vignette describes the two latent variable models implemented in
`lblvm`: the **bounded-continuous response model** and the **discretized
response model (DRM)**. Both models use the same logit-beta framework,
but differ in how the response is observed.

![\*\*Figure 1\*\* Bounded-continuous response model in Panel (a) and
discretized response model in Panel (b). The black solid lines in the
two panels are identical and are the \$\mu\$-curves. The grey curves
represent the probability densities of continuous item responses in
Panel (a) and those of internal judgment continuum in Panel (b). The
grey line segments, located exactly on the \$\mu\$-curve, indicate the
means of the densities. The relative width of the color bands in Panel
(b) can be estimated in the discretized response
model.](Models_files/figure-html/unnamed-chunk-2-1.png)

**Figure 1** Bounded-continuous response model in Panel (a) and
discretized response model in Panel (b). The black solid lines in the
two panels are identical and are the $`\mu`$-curves. The grey curves
represent the probability densities of continuous item responses in
Panel (a) and those of internal judgment continuum in Panel (b). The
grey line segments, located exactly on the $`\mu`$-curve, indicate the
means of the densities. The relative width of the color bands in Panel
(b) can be estimated in the discretized response model.

## Bounded-continuous response model

The bounded-continuous response model is designed for responses that
take values between 0 and 1. Let $`x_{ij}`$ denote the response of
person $`i`$ to item $`j`$, and let $`\theta_i`$ denote the person’s
latent trait. The expected response is modeled as

``` math
\operatorname{logit}(\mu_{ij})
=
a_j\theta_i+c_j,
```

where

``` math
\mu_{ij}
=
E(x_{ij}\mid\theta_i)
=
\frac{\exp(a_j\theta_i+c_j)}
{1+\exp(a_j\theta_i+c_j)}.
```

Here, $`a_j`$ is the item slope and $`c_j`$ is the item intercept. The
slope determines how rapidly the expected response changes as the latent
trait changes, whereas the intercept determines the location of the
expected response curve. The latent trait can be multidimensional, and
in such cases, $`a`$ and $`\theta`$ become $`D`$-dimensional vectors and
$`a_j\theta_i`$ in the above equation is substituted with an inner
product
$`\mathbf{a}_j\cdot\mathbf{\theta}_i=\sum_{d=1}^{D}a_{jd}\theta_{id}`$.

The response is then assumed to follow a beta distribution,

``` math
x_{ij}\mid\theta_i
\sim
\operatorname{Beta}
\left(
\mu_{ij}\nu_j,
(1-\mu_{ij})\nu_j
\right),
```

where $`\nu_j>0`$ is an item-specific precision parameter. Under this
parameterization,

``` math
\operatorname{Var}(x_{ij}\mid\theta_i)
=
\frac{\mu_{ij}(1-\mu_{ij})}{\nu_j+1}.
```

The precision parameter controls how concentrated the responses are
around their expected value. Larger values of $`\nu_j`$ correspond to
greater precision and, therefore, less residual variation.

The bounded-continuous model can therefore be understood in two parts.
The logit equation describes the relationship between the latent trait
and the expected response, while the beta distribution describes how
individual responses vary around that expectation.

## Discretized response model

The **discretized response model (DRM)** extends the bounded-continuous
model to ordered categorical responses. The key assumption is that an
ordered categorical response arises by discretizing an underlying
bounded-continuous judgment.

For example, consider an item with five response categories:

``` text
Strongly disagree | Disagree | Neutral | Agree | Strongly agree
```

DRM assumes that these categories correspond to five consecutive
intervals on an underlying continuum from 0 to 1. Let the boundaries of
these intervals be

``` math
0=T_0<T_1<T_2<T_3<T_4<T_5=1.
```

The underlying continuous response follows the same beta distribution as
the bounded-continuous model. Therefore, for category $`z`$, the
probability of observing that category is the probability that the
underlying response falls within its corresponding interval:

``` math
P(U_{ij}=z\mid\theta_i)
=
\int_{T_{z-1}}^{T_z}
\operatorname{Beta}
\left(x \mid
\mu_{ij}\nu_j,
(1-\mu_{ij})\nu_j
\right)\,dx,
```

where
$`\operatorname{Beta}\left(x \mid \mu_{ij}\nu_j,(1-\mu_{ij})\nu_j\right)`$
is the beta density implied by the bounded-continuous model.

For example, the probability of the fourth category is

``` math
P(U_{ij}=4\mid\theta_i)
=
\int_{T_3}^{T_4}
\operatorname{Beta}
\left(x \mid
\mu_{ij}\nu_j,
(1-\mu_{ij})\nu_j
\right)\,dx.
```

Thus, the probability of a category is determined by both the underlying
judgment distribution and the interval assigned to that category.

## Relationship between the two models

The two models share the same underlying logit-beta response process.
Their difference is the observation mechanism.

For the bounded-continuous response model, the underlying response is
observed directly:

``` math
\theta_i
\longrightarrow
\mu_{ij}
\longrightarrow
x_{ij}.
```

For DRM, the same underlying response is discretized into ordered
categories:

``` math
\theta_i
\longrightarrow
\mu_{ij}
\longrightarrow
x_{ij}
\longrightarrow
U_{ij}.
```

The first two stages are identical. The difference is that DRM
introduces thresholds that divide the bounded continuum into intervals,
and the observed category indicates which interval contains the
underlying judgment.

Consequently, DRM allows the latent-trait effects and the response
categorization process to be represented separately: the latent trait
determines the location and precision of the underlying judgment
distribution, while the thresholds determine how that judgment is
translated into an ordered response.
