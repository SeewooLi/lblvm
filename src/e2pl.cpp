#include <RcppArmadillo.h>
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
// [[Rcpp::export]]
List cont_L1L2_cpp(const arma::vec& item,
                   const arma::mat& grid,
                   const arma::vec& data,
                   const arma::mat& Pk,
                   bool calculate_m = false) {

  const int N = data.n_elem;
  const int M = grid.n_rows;
  const int D = grid.n_cols;

  const double P_EPS = 1e-12;


  const arma::vec a = item.head(D);
  const double c = item[D];

  const double xi = item[D + 1];
  const double nu = std::exp(xi);

  // nu does not depend on m
  const double digamma_nu  = R::digamma(nu);
  const double trigamma_nu = R::trigamma(nu);


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

  arma::vec L1(D + 2, fill::zeros);
  arma::mat H(D + 2, D + 2, fill::zeros);


  for (int m = 0; m < M; ++m) {

    // eta = theta'a + c
    const arma::rowvec theta = grid.row(m);
    const double eta = arma::dot(theta, a) + c;

    // Stable logistic transformation
    double mu;

    if (eta >= 0.0) {
      const double e = std::exp(-eta);
      mu = 1.0 / (1.0 + e);
    } else {
      const double e = std::exp(eta);
      mu = e / (1.0 + e);
    }

    const double q = mu * (1.0 - mu);

    const double alpha = nu * mu;
    const double beta  = nu * (1.0 - mu);

    const double da = R::digamma(alpha);
    const double db = R::digamma(beta);

    const double ta = R::trigamma(alpha);
    const double tb = R::trigamma(beta);


    // Derivatives with respect to eta

    const double S =
      s1[m]
    - s2[m]
    - f[m] * (da - db);

    const double L_eta =
    nu * q * S;

    const double L_eta_eta =
      nu * q * (1.0 - 2.0 * mu) * S
    - (nu * q) * (nu * q)
      * f[m] * (ta + tb);


    const double L_xi =
    nu * (
        f[m] * digamma_nu
    + mu * (s1[m] - f[m] * da)
      + (1.0 - mu) * (s2[m] - f[m] * db)
    );

    const double L_eta_xi =
      nu * q * (
          S
          - f[m] * (
              alpha * ta
    - beta * tb
          )
      );

    const double L_xi_xi =
      L_xi
      + nu * nu * f[m] * trigamma_nu
    - f[m] * (
        alpha * alpha * ta
    + beta * beta * tb
    );


    // Gradient

    for (int j = 0; j < D; ++j)
      L1[j] += theta[j] * L_eta;

    L1[D] += L_eta;

    L1[D + 1] += L_xi;


    // Hessian

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


  // Missing information matrix

  if (calculate_m) {

    arma::mat IMm;
    IMm.zeros(D + 2, D + 2);

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
      arma::vec score_mean(D + 2, fill::zeros);

      // E[S S' | x]
      arma::mat score_second(D + 2, D + 2, fill::zeros);

      for (int m = 0; m < M; ++m) {

        const double post = Pk(i, m);

        if (post == 0.0)
          continue;

        const arma::rowvec theta = grid.row(m);

        const double eta =
          arma::dot(theta, a) + c;

        // Stable logistic transformation
        double mu;

        if (eta >= 0.0) {
          const double e = std::exp(-eta);
          mu = 1.0 / (1.0 + e);
        } else {
          const double e = std::exp(eta);
          mu = e / (1.0 + e);
        }

        const double alpha = nu * mu;
        const double beta  = nu * (1.0 - mu);

        const double da = R::digamma(alpha);
        const double db = R::digamma(beta);

        const double w1 = log_x  - da;
        const double w2 = log_1x - db;

        const double temp = w1 - w2;

        // Complete-data score conditional on latent m
        arma::vec score(D + 2);

        for (int d = 0; d < D; ++d)
          score[d] = theta[d] * temp;

        score[D] = temp;

        // Score with respect to xi = log(nu)
        score[D + 1] =
          nu * (
              digamma_nu
              + mu * w1
        + (1.0 - mu) * w2
          );

        score_mean += post * score;
        score_second += post * (score * score.t());
      }

      // Var(S | x)
      IMm += score_second
      - score_mean * score_mean.t();
    }
    return List::create(
      _["gradient"] = L1,
      _["hessian"]  = H,
      _["IMm"]      = IMm
    );

  } else {

    return List::create(
      _["gradient"] = L1,
      _["hessian"]  = H
    );
  }
}
