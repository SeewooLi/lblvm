#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
List Estep_Cont_cpp(const NumericMatrix& data,
                    const NumericMatrix& item,
                    const NumericMatrix& grid,
                    const NumericVector& prior) {

  const int N = data.nrow();
  const int I = data.ncol();
  const int M = grid.nrow();
  const int D = grid.ncol();

  const double P_EPS = 1e-12;

  // Precompute sigmoid for every item × quadrature point

  NumericMatrix p(I, M);

  for (int i = 0; i < I; i++) {

    for (int m = 0; m < M; m++) {

      double eta = item(i, D);  // c

      for (int d = 0; d < D; d++)
        eta += grid(m, d) * item(i, d);

      // expit(eta)
      p(i, m) = 1.0 / (1.0 + std::exp(-eta));

      p(i, m) = std::min(
        std::max(p(i, m), P_EPS),
        1.0 - P_EPS
      );
    }
  }


  // Person-quadrature-point log likelihood

  NumericMatrix logL(N, M);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int n = 0; n < N; n++) {

    for (int m = 0; m < M; m++) {

      double ll = 0.0;

      for (int i = 0; i < I; i++) {

        double x = data(n, i);

        x = std::min(
          std::max(x, P_EPS),
          1.0 - P_EPS
        );

        // Missing response
        if (NumericVector::is_na(x))
          continue;

        double pp = p(i, m);
        double nu = std::exp(item(i, D + 1));

        double alpha = pp * nu;
        double beta  = (1.0 - pp) * nu;

        ll += std::lgamma(nu)
          - std::lgamma(alpha)
          - std::lgamma(beta)
          + (alpha - 1.0) * std::log(x)
          + (beta - 1.0) * std::log1p(-x);
      }

      logL(n, m) = ll;
    }
  }


  // Posterior probabilities

  NumericMatrix posterior(N, M);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int n = 0; n < N; n++) {

    double denom = 0.0;

    for (int m = 0; m < M; m++) {

      if (prior[m] > 0.0) {

        posterior(n, m) = std::exp(logL(n, m) + std::log(prior[m]));

        denom += posterior(n, m);

      } else {

        posterior(n, m) = 0.0;
      }
    }

    // Normalize
    for (int m = 0; m < M; m++)
      posterior(n, m) /= denom;
  }


  // Expected frequency

  NumericVector freq(M);
  double total_freq = 0.0;

  for (int m = 0; m < M; m++) {

    double sum = 0.0;

    for (int n = 0; n < N; n++)
      sum += posterior(n, m);

    freq[m] = sum;
    total_freq += freq[m];
  }

  NumericVector Ak(M);

  for (int m = 0; m < M; m++)
    Ak[m] = freq[m] / total_freq;


  // log likelihood

  double total_logL = 0.0;

  for (int n = 0; n < N; n++) {

    for (int m = 0; m < M; m++) {

      double r = posterior(n, m);

      if (r > 0.0) {

        total_logL += r * (
          logL(n, m)
        + std::log(prior[m])
        - std::log(r)
        );
      }
    }
  }


  return List::create(
    _["posterior"] = posterior,
    _["freq"]      = freq,
    _["grid"]      = grid,
    _["prior"]     = prior,
    _["logL"]      = total_logL,
    _["Ak"]        = Ak
  );
}
/*** R
# data
data <- matrix(runif(400), ncol = 4)

# parameter
item <- cbind(1,1,0,rep(exp(2),4))

# grid
x <- seq(-4, 4, length.out=41)
grid_list <- replicate(2, x, simplify = FALSE)
grid <- as.matrix(do.call(expand.grid, grid_list))

# prior
prior <- mvtnorm::dmvnorm(grid,
                          mean = c(0,0),
                          sigma = diag(c(1,1))
)
prior <- prior/sum(prior)

aaa <- Estep_Cont_cpp(data, item, grid, prior)

x <- sort(unique(grid[, "Var1"]))
y <- sort(unique(grid[, "Var2"]))
z <- matrix(aaa$posterior[1,], nrow = length(x), ncol = length(y))

#persp(x, y, z,
#      theta = 30, phi = 30,          # Viewing angles (azimuth and colatitude)
#      expand = 0.5,                  # Scale factor for height axis
#      col = "lightblue",             # Surface color
#      shade = 0.5,                   # Add shading effect
#      xlab = "Var1", ylab = "Var2", zlab = "Posterior")
*/


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
List cont_L1L2_cpp(const arma::vec& item,
                   const arma::mat& grid,
                   const arma::vec& data,
                   const arma::mat& Pk,
                   bool calculate_m = false) {

  const int N = data.n_elem;
  const int M = grid.n_rows;
  const int D = grid.n_cols;
  const int P = D + 2;

  const double P_EPS = 1e-12;

  const arma::vec a = item.head(D);
  const double c = item[D];

  const double xi = item[D + 1];
  const double nu = std::exp(xi);

  const double digamma_nu  = R::digamma(nu);
  const double trigamma_nu = R::trigamma(nu);


  arma::vec mu(M);
  arma::vec q(M);

  arma::vec alpha(M);
  arma::vec beta(M);

  arma::vec da(M);
  arma::vec db(M);

  arma::vec ta(M);
  arma::vec tb(M);

  for (int m = 0; m < M; ++m) {

    const double eta =
      arma::dot(grid.row(m), a) + c;

    // Stable logistic transformation
    if (eta >= 0.0) {

      const double e = std::exp(-eta);
      mu[m] = 1.0 / (1.0 + e);

    } else {

      const double e = std::exp(eta);
      mu[m] = e / (1.0 + e);
    }

    q[m] = mu[m] * (1.0 - mu[m]);

    alpha[m] = nu * mu[m];
    beta[m]  = nu * (1.0 - mu[m]);

    da[m] = R::digamma(alpha[m]);
    db[m] = R::digamma(beta[m]);

    ta[m] = R::trigamma(alpha[m]);
    tb[m] = R::trigamma(beta[m]);
  }


  // Sufficient statistics from the E-step

  arma::vec f(M, fill::zeros);
  arma::vec s1(M, fill::zeros);
  arma::vec s2(M, fill::zeros);

  for (int n = 0; n < N; ++n) {

    if (!std::isfinite(data[n]))
      continue;

    double x = data[n];

    x = std::min(
      std::max(x, P_EPS),
      1.0 - P_EPS
    );

    const double log_x  = std::log(x);
    const double log_1x = std::log1p(-x);

    for (int m = 0; m < M; ++m) {

      const double w = Pk(n, m);

      f[m]  += w;
      s1[m] += w * log_x;
      s2[m] += w * log_1x;
    }
  }


  // Gradient and Hessian

  arma::vec L1(P, fill::zeros);
  arma::mat H(P, P, fill::zeros);

  for (int m = 0; m < M; ++m) {

    const arma::rowvec theta = grid.row(m);

    const double S =
      s1[m]
    - s2[m]
    - f[m] * (da[m] - db[m]);

    const double L_eta =
    nu * q[m] * S;

    const double L_eta_eta =
      nu * q[m] * (1.0 - 2.0 * mu[m]) * S
    - (nu * q[m]) * (nu * q[m])
      * f[m] * (ta[m] + tb[m]);

      const double L_xi =
      nu * (
          f[m] * digamma_nu
      + mu[m] * (s1[m] - f[m] * da[m])
        + (1.0 - mu[m]) * (s2[m] - f[m] * db[m])
      );

      const double L_eta_xi =
        nu * q[m] * (
            S
            - f[m] * (
                alpha[m] * ta[m]
      - beta[m] * tb[m]
            )
        );

      const double L_xi_xi =
        L_xi
        + nu * nu * f[m] * trigamma_nu
      - f[m] * (
          alpha[m] * alpha[m] * ta[m]
      + beta[m] * beta[m] * tb[m]
      );

      for (int j = 0; j < D; ++j)
        L1[j] += theta[j] * L_eta;

      L1[D] += L_eta;
      L1[D + 1] += L_xi;


      for (int j = 0; j < D; ++j) {

        for (int k = 0; k < D; ++k)
          H(j, k) +=
            L_eta_eta * theta[j] * theta[k];

        H(j, D) +=
          L_eta_eta * theta[j];

        H(D, j) +=
          L_eta_eta * theta[j];
      }

      H(D, D) += L_eta_eta;

      for (int j = 0; j < D; ++j) {

        H(j, D + 1) +=
          L_eta_xi * theta[j];

        H(D + 1, j) +=
          L_eta_xi * theta[j];
      }

      H(D, D + 1) += L_eta_xi;
      H(D + 1, D) += L_eta_xi;

      H(D + 1, D + 1) += L_xi_xi;
  }


  if (!calculate_m) {

    return List::create(
      _["gradient"] = L1,
      _["hessian"]  = H
    );
  }


  // Missing information matrix

  arma::mat IMm(P, P, fill::zeros);

#ifdef _OPENMP

  const int n_threads = omp_get_max_threads();

  std::vector<arma::mat> IMm_thread(
      n_threads,
      arma::mat(P, P, fill::zeros)
  );

#pragma omp parallel
{

  const int tid = omp_get_thread_num();

  arma::mat& IMm_local = IMm_thread[tid];

#pragma omp for schedule(static)

  for (int i = 0; i < N; ++i) {

    if (!std::isfinite(data[i]))
      continue;

    double x = data[i];

    x = std::min(
      std::max(x, P_EPS),
      1.0 - P_EPS
    );

    const double log_x  = std::log(x);
    const double log_1x = std::log1p(-x);


    // E[S | x]

    arma::vec score_mean(P, fill::zeros);

    // E[S S' | x]

    arma::mat score_second(P, P, fill::zeros);


    for (int m = 0; m < M; ++m) {

      const double post = Pk(i, m);

      if (post == 0.0)
        continue;

      const arma::rowvec theta = grid.row(m);


      const double w1 =
        log_x - da[m];

      const double w2 =
        log_1x - db[m];

      const double temp =
        w1 - w2;

      const double score_xi =
        nu * (
            digamma_nu
            + mu[m] * w1
      + (1.0 - mu[m]) * w2
        );


      const double wt =
        post * temp;

      for (int d = 0; d < D; ++d)
        score_mean[d] +=
          wt * theta[d];

      score_mean[D] += wt;

      score_mean[D + 1] +=
        post * score_xi;


      const double wtt =
        post * temp * temp;

      const double wtx =
        post * temp * score_xi;

      const double wxx =
        post * score_xi * score_xi;


      for (int d = 0; d < D; ++d) {

        for (int k = 0; k < D; ++k) {

          score_second(d, k) +=
            wtt * theta[d] * theta[k];
        }
      }


      for (int d = 0; d < D; ++d) {

        score_second(d, D) +=
          wtt * theta[d];

        score_second(D, d) +=
          wtt * theta[d];
      }


      for (int d = 0; d < D; ++d) {

        score_second(d, D + 1) +=
          wtx * theta[d];

        score_second(D + 1, d) +=
          wtx * theta[d];
      }


      score_second(D, D) +=
        wtt;

      score_second(D, D + 1) +=
        wtx;

      score_second(D + 1, D) +=
        wtx;

      score_second(D + 1, D + 1) +=
        wxx;
    }


    IMm_local +=
      score_second
      - score_mean * score_mean.t();
  }
}


// Reduce thread-local matrices

for (int t = 0; t < n_threads; ++t)
  IMm += IMm_thread[t];

#else

// Serial fallback

for (int i = 0; i < N; ++i) {

  if (!std::isfinite(data[i]))
    continue;

  double x = data[i];

  x = std::min(
    std::max(x, P_EPS),
    1.0 - P_EPS
  );

  const double log_x  = std::log(x);
  const double log_1x = std::log1p(-x);


  arma::vec score_mean(P, fill::zeros);
  arma::mat score_second(P, P, fill::zeros);


  for (int m = 0; m < M; ++m) {

    const double post = Pk(i, m);

    if (post == 0.0)
      continue;

    const arma::rowvec theta = grid.row(m);

    const double w1 =
      log_x - da[m];

    const double w2 =
      log_1x - db[m];

    const double temp =
      w1 - w2;

    const double score_xi =
      nu * (
          digamma_nu
          + mu[m] * w1
    + (1.0 - mu[m]) * w2
      );


    // ------------------------------------------------------------
    // E[S | x]
    // ------------------------------------------------------------

    const double wt =
      post * temp;

    for (int d = 0; d < D; ++d)
      score_mean[d] +=
        wt * theta[d];

    score_mean[D] += wt;

    score_mean[D + 1] +=
      post * score_xi;


    // ------------------------------------------------------------
    // E[S S' | x]
    // ------------------------------------------------------------

    const double wtt =
      post * temp * temp;

    const double wtx =
      post * temp * score_xi;

    const double wxx =
      post * score_xi * score_xi;


    for (int d = 0; d < D; ++d) {

      for (int k = 0; k < D; ++k) {

        score_second(d, k) +=
          wtt * theta[d] * theta[k];
      }
    }


    for (int d = 0; d < D; ++d) {

      score_second(d, D) +=
        wtt * theta[d];

      score_second(D, d) +=
        wtt * theta[d];

      score_second(d, D + 1) +=
        wtx * theta[d];

      score_second(D + 1, d) +=
        wtx * theta[d];
    }


    score_second(D, D) +=
      wtt;

    score_second(D, D + 1) +=
      wtx;

    score_second(D + 1, D) +=
      wtx;

    score_second(D + 1, D + 1) +=
      wxx;
  }


  IMm +=
    score_second
    - score_mean * score_mean.t();
}

#endif


// ================================================================
// Return
// ================================================================

return List::create(
  _["gradient"] = L1,
  _["hessian"]  = H,
  _["IMm"]      = IMm
);
}
