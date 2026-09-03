// Koop and Korobilis (2014), "A new index of financial conditions".
//
// A time varying parameter FAVAR estimated without simulation. The loadings of
// the measurement equation and the coefficients of the transition equation both
// drift as random walks, but neither the variance of that drift nor the two
// error covariances is given a prior and drawn: they are replaced by forgetting
// factors, which is what makes the estimator a pair of Kalman passes rather than
// a Gibbs sampler. See the package vignette for what that costs.
//
// The state is z_t = [y_t', f_t']', the observed macro block on top of the
// latent factors, and the observation is [y_t', x_t']'. The macro block is
// measured without error -- the leading n_y diagonal entries of the measurement
// covariance are zero and the leading n_y rows of the loading matrix are
// [I 0] -- so the filter reproduces y_t in the state exactly and only the
// factors are inferred.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

namespace
{

// Solves F X = rhs' for a small symmetric F, falling back on the pseudo-inverse
// when F is rank deficient.
//
// It regularly is, and structurally rather than by accident: the macro block is
// measured without error, so once a period has been seen its state has no
// conditional variance left, and the lag blocks of the predicted state
// covariance are exactly singular from the second period on. The pseudo-inverse
// is the right answer there rather than a rescue -- a direction that is known
// exactly needs no correction -- so the approximate solution Armadillo would
// otherwise warn about is suppressed with no_approx and taken deliberately here.
arma::mat right_solve(const arma::mat &rhs, const arma::mat &f)
{
  arma::mat out;
  if (arma::solve(out, f, rhs.t(),
                  arma::solve_opts::likely_sympd + arma::solve_opts::no_approx))
  {
    return out.t();
  }
  return rhs * arma::pinv(f);
}

// log density of a multivariate normal with covariance `f`, evaluated at the
// prediction error `e`. Used for the one step ahead predictive likelihood, which
// is what a dynamic model averaging layer weights the model by. Returns NA when
// the covariance is not positive definite rather than throwing: a single
// unusable period should not lose the whole path.
double log_density(const arma::vec &e, const arma::mat &f)
{
  double sign = 0.0;
  double log_det = 0.0;
  if (!arma::log_det(log_det, sign, f) || sign <= 0.0)
  {
    return NA_REAL;
  }
  arma::vec sol;
  if (!arma::solve(sol, f, e,
                   arma::solve_opts::likely_sympd + arma::solve_opts::no_approx))
  {
    return NA_REAL;
  }
  return -0.5 * (static_cast<double>(e.n_elem) * std::log(2.0 * M_PI) + log_det +
                 arma::dot(e, sol));
}

// The companion matrix of the transition, built from one period's coefficients.
// `beta` is vec(B) with B of n_z x (n_z p + n_const); the intercept column, when
// there is one, is the last and is returned separately in `intercept`.
arma::mat companion(const arma::vec &beta, const int n_z, const int p,
                    const int n_const, arma::vec &intercept)
{
  const int n_reg = n_z * p + n_const;
  const arma::mat b = arma::reshape(beta, n_z, n_reg);

  const int n_s = n_z * p;
  arma::mat comp(n_s, n_s, arma::fill::zeros);
  comp.rows(0, n_z - 1) = b.cols(0, n_s - 1);
  if (p > 1)
  {
    comp.submat(n_z, 0, n_s - 1, n_s - n_z - 1) =
      arma::eye<arma::mat>(n_s - n_z, n_s - n_z);
  }

  intercept.zeros(n_s);
  if (n_const > 0)
  {
    intercept.head(n_z) = b.col(n_reg - 1);
  }
  return comp;
}

// The regressors of the transition equation, one column per period:
// [z_{t-1}' .. z_{t-p}' (1)]'. Periods before the sample are zero, and the
// filter does not update the coefficients until period p anyway.
arma::mat transition_regressors(const arma::mat &z, const int p, const int n_const)
{
  const arma::uword tt = z.n_rows;
  const arma::uword n_z = z.n_cols;
  arma::mat out(n_z * p + n_const, tt, arma::fill::zeros);
  for (int j = 1; j <= p; j++)
  {
    for (arma::uword t = static_cast<arma::uword>(j); t < tt; t++)
    {
      out.submat((j - 1) * n_z, t, j * n_z - 1, t) = z.row(t - j).t();
    }
  }
  if (n_const > 0)
  {
    out.row(n_z * p + n_const - 1).ones();
  }
  return out;
}

// Principal components of an already normalised matrix, scaled the way Koop and
// Korobilis scale them: the loadings are sqrt(n_x) times the eigenvectors, and
// the factors are x L / n_x. The scale is arbitrary in a factor model but it is
// the starting value the filter runs from, so it is worth matching.
arma::mat principal_components(const arma::mat &x, const int n_fac)
{
  arma::mat coeff;
  arma::vec eigval;
  const arma::mat xx = arma::cov(x);
  if (!arma::eig_sym(eigval, coeff, xx))
  {
    Rcpp::stop("The principal component decomposition of 'x' failed.");
  }
  // eig_sym returns ascending eigenvalues; take the largest n_fac.
  const arma::uword n_x = x.n_cols;
  arma::mat load(n_x, n_fac);
  for (int j = 0; j < n_fac; j++)
  {
    load.col(j) = coeff.col(n_x - 1 - static_cast<arma::uword>(j));
  }
  load *= std::sqrt(static_cast<double>(n_x));
  return x * load / static_cast<double>(n_x);
}

} // namespace

//' Koop and Korobilis (2014) TVP-FAVAR
//'
//' The estimator behind \code{\link{add_posterior_coefficients.fcimodel}}. Not
//' exported: it takes the model object apart itself and assumes
//' \code{\link{add_priors.fcimodel}} has filled it in.
//'
//' @param object an object of class \code{"fcimodel"}.
//'
//' @noRd
// [[Rcpp::export(.kk_favar)]]
Rcpp::List kk_favar(Rcpp::List object)
{
  const Rcpp::List data = object["data"];
  const Rcpp::List model = object["model"];
  const Rcpp::List priors = object["priors"];

  const arma::mat y = Rcpp::as<arma::mat>(data["y"]); // tt x n_y
  const arma::mat x = Rcpp::as<arma::mat>(data["x"]); // tt x n_x

  const int n_fac = Rcpp::as<int>(model["n"]);
  const int p = Rcpp::as<int>(model["p"]);
  const int n_const = Rcpp::as<bool>(model["constant"]) ? 1 : 0;
  const int passes = Rcpp::as<int>(model["iterations"]);

  const Rcpp::List decay = model["decay"];
  const double l_u = Rcpp::as<double>(decay["u"]);           // EWMA on V
  const double l_v = Rcpp::as<double>(decay["v"]);           // EWMA on Q
  const double l_lambda = Rcpp::as<double>(decay["lambda"]); // discount on the loading state
  const double l_beta = Rcpp::as<double>(decay["beta"]);     // discount on the transition state

  const arma::uword tt = y.n_rows;
  const arma::uword n_y = y.n_cols;
  const arma::uword n_x = x.n_cols;
  const arma::uword n_z = n_y + static_cast<arma::uword>(n_fac);
  const arma::uword n_s = n_z * static_cast<arma::uword>(p);
  const arma::uword n_reg = n_z * static_cast<arma::uword>(p) + static_cast<arma::uword>(n_const);
  const arma::uword n_beta = n_z * n_reg;
  const arma::uword q = n_y + n_x; // width of the observation vector

  // Priors -- initial conditions only. Nothing here is a prior on a variance of
  // the drift or of an error: those are the four forgetting factors above.
  const Rcpp::List prior_lambda = priors["lambda"];
  const Rcpp::List prior_beta = priors["beta"];
  const Rcpp::List prior_factors = priors["factors"];
  const arma::vec lambda_mu = Rcpp::as<arma::vec>(prior_lambda["mu"]);   // n_z
  const arma::mat lambda_v = Rcpp::as<arma::mat>(prior_lambda["v"]);     // n_z x n_z
  const arma::vec beta_mu = Rcpp::as<arma::vec>(prior_beta["mu"]);       // n_beta
  const arma::mat beta_v = Rcpp::as<arma::mat>(prior_beta["v"]);         // n_beta x n_beta
  const arma::vec factor_mu = Rcpp::as<arma::vec>(prior_factors["mu"]);  // n_s
  const arma::mat factor_v = Rcpp::as<arma::mat>(prior_factors["v"]);    // n_s x n_s
  const arma::vec u_init = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(priors["u"])["v"]); // n_x
  const arma::vec v_init = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(priors["v"])["v"]); // n_z

  // Starting factor path. add_priors() puts the principal components here; the
  // C++ recomputes them only if it was handed nothing.
  arma::mat factors;
  if (data.containsElementNamed("factors") &&
      !Rf_isNull(data["factors"]))
  {
    factors = Rcpp::as<arma::mat>(data["factors"]);
  }
  else
  {
    factors = principal_components(x, n_fac);
  }

  // Paths carried between the two steps and returned.
  arma::mat lambda(tt, n_x * n_z, arma::fill::zeros); // row-major within a period
  arma::mat beta(tt, n_beta, arma::fill::zeros);
  arma::mat u_sigma(tt, n_x, arma::fill::zeros); // diagonal of V_t
  arma::mat v_sigma(tt, n_z, arma::fill::zeros); // diagonal of Q_t
  arma::vec loglik(tt, arma::fill::value(NA_REAL));

  for (int pass = 0; pass < passes; pass++)
  {
    Rcpp::checkUserInterrupt();

    // The state as it currently stands: the observed macro block beside the
    // factors of the previous pass.
    arma::mat z(tt, n_z);
    z.cols(0, n_y - 1) = y;
    z.cols(n_y, n_z - 1) = factors;

    const arma::mat z_lag = transition_regressors(z, p, n_const);

    // ---------------------------------------------------------------------
    // Step one: the parameters, given the factors. Koop and Korobilis'
    // KFS_parameters.
    // ---------------------------------------------------------------------

    // The loadings, row by row. Row i of x has its own n_z dimensional state
    // against the design z_t, and given a diagonal V the rows are independent.
    arma::mat lambda_upd(tt, n_x * n_z, arma::fill::zeros);
    {
      arma::mat state(n_x, n_z);
      for (arma::uword i = 0; i < n_x; i++) { state.row(i) = lambda_mu.t(); }
      arma::cube s_cov(n_z, n_z, n_x);
      for (arma::uword i = 0; i < n_x; i++) { s_cov.slice(i) = lambda_v; }

      arma::vec v_diag = u_init;

      for (arma::uword t = 0; t < tt; t++)
      {
        const arma::vec z_t = z.row(t).t();
        const arma::vec x_t = x.row(t).t();

        // The error covariance is discounted from the prediction errors before
        // it enters the gain, which is the order the paper's code uses.
        const arma::vec e = x_t - state * z_t;

        for (arma::uword i = 0; i < n_x; i++)
        {
          // Predicted state covariance: the random walk innovation is not a
          // parameter here, it is whatever the forgetting factor implies.
          const arma::mat r = (t == 0) ? s_cov.slice(i) : s_cov.slice(i) / l_lambda;

          // A series that has not started yet has nothing to update its loading
          // or its error variance with. The loading still drifts, though, so the
          // prediction stands as this period's estimate and the covariance keeps
          // the period's discounting -- which is what the backward pass, where
          // every period contributes its own factor of kappa, assumes has
          // happened.
          if (!std::isfinite(x_t(i)))
          {
            s_cov.slice(i) = r;
            continue;
          }
          if (t > 0)
          {
            v_diag(i) = l_u * v_diag(i) + (1.0 - l_u) * e(i) * e(i);
          }

          const arma::vec rz = r * z_t;
          const double f = v_diag(i) + arma::dot(z_t, rz);
          if (f > 0.0)
          {
            const arma::vec k = rz / f;
            state.row(i) += (k * e(i)).t();
            s_cov.slice(i) = r - k * rz.t();
          }
          else
          {
            s_cov.slice(i) = r;
          }
        }
        u_sigma.row(t) = v_diag.t();
        lambda_upd.row(t) = arma::vectorise(state, 1);
      }

      // The backward pass. For a random walk whose predicted covariance is
      // S_{t-1} / l, the smoother gain S_t R_{t+1}^{-1} is exactly l times the
      // identity, so the Rauch-Tung-Striebel recursion collapses to backward
      // exponential smoothing. No covariance path has to be stored for it.
      lambda = lambda_upd;
      for (arma::uword t = tt - 1; t-- > 0;)
      {
        lambda.row(t) = (1.0 - l_lambda) * lambda_upd.row(t) + l_lambda * lambda.row(t + 1);
      }
    }

    // The transition coefficients, as one state of n_beta elements against the
    // SUR design kron(z_lag_t', I).
    arma::mat beta_upd(tt, n_beta, arma::fill::zeros);
    {
      arma::vec state = beta_mu;
      arma::mat s_cov = beta_v;
      arma::vec q_diag = v_init;
      const arma::mat eye_z = arma::eye<arma::mat>(n_z, n_z);

      for (arma::uword t = 0; t < tt; t++)
      {
        // Nothing to learn from before the first complete set of lags.
        if (t >= static_cast<arma::uword>(p))
        {
          const arma::mat design = arma::kron(z_lag.col(t).t(), eye_z); // n_z x n_beta
          const arma::vec e = z.row(t).t() - design * state;

          q_diag = l_v * q_diag + (1.0 - l_v) * arma::square(e);

          const arma::mat r = s_cov / l_beta;
          const arma::mat rz = r * design.t(); // n_beta x n_z
          const arma::mat f = design * rz + arma::diagmat(q_diag);
          const arma::mat k = right_solve(rz, f); // n_beta x n_z
          state += k * e;
          s_cov = r - k * rz.t();
        }
        beta_upd.row(t) = state.t();
        v_sigma.row(t) = q_diag.t();
      }

      beta = beta_upd;
      for (arma::uword t = tt - 1; t-- > 0;)
      {
        beta.row(t) = (1.0 - l_beta) * beta_upd.row(t) + l_beta * beta.row(t + 1);
      }
    }

    // ---------------------------------------------------------------------
    // Step two: the factors, given the parameters. Koop and Korobilis'
    // KFS_factors, on the companion form of the transition.
    // ---------------------------------------------------------------------

    arma::mat s_pred(n_s, tt, arma::fill::zeros);
    arma::mat s_upd(n_s, tt, arma::fill::zeros);
    arma::cube r_cov(n_s, n_s, tt, arma::fill::zeros);
    arma::cube s_cov(n_s, n_s, tt, arma::fill::zeros);
    arma::cube comp(n_s, n_s, tt, arma::fill::zeros);

    // The measurement matrix. Its leading n_y rows are [I 0] in every period and
    // are never estimated: the macro block of the state is the observed data.
    arma::mat load(q, n_z, arma::fill::zeros);
    load.submat(0, 0, n_y - 1, n_y - 1) = arma::eye<arma::mat>(n_y, n_y);

    for (arma::uword t = 0; t < tt; t++)
    {
      arma::vec intercept(n_s, arma::fill::zeros);
      comp.slice(t) = companion(beta.row(t).t(), n_z, p, n_const, intercept);

      if (t == 0)
      {
        s_pred.col(t) = factor_mu;
        r_cov.slice(t) = factor_v;
      }
      else
      {
        s_pred.col(t) = comp.slice(t - 1) * s_upd.col(t - 1) + intercept;
        arma::mat state_var(n_s, n_s, arma::fill::zeros);
        state_var.submat(0, 0, n_z - 1, n_z - 1) = arma::diagmat(v_sigma.row(t).t());
        r_cov.slice(t) =
          comp.slice(t - 1) * s_cov.slice(t - 1) * comp.slice(t - 1).t() + state_var;
      }

      load.submat(n_y, 0, q - 1, n_z - 1) =
        arma::reshape(lambda.row(t), n_z, n_x).t();

      // The measurement matrix over the whole companion state. Only the current
      // period is observed, so the lag blocks are zero -- but the update still
      // has to run over the full state, because that is what propagates the
      // information back into the lags the next period predicts from.
      arma::mat design(q, n_s, arma::fill::zeros);
      design.cols(0, n_z - 1) = load;

      arma::vec obs(q);
      obs.head(n_y) = y.row(t).t();
      obs.tail(n_x) = x.row(t).t();

      // The macro block of the measurement covariance is zero: y is observed.
      arma::mat meas_var(q, q, arma::fill::zeros);
      meas_var.submat(n_y, n_y, q - 1, q - 1) = arma::diagmat(u_sigma.row(t).t());

      // Rows without an observation are dropped from this period's measurement
      // rather than filled in. The macro block is always present.
      arma::uvec seen(q);
      arma::uword n_seen = 0;
      for (arma::uword j = 0; j < q; j++)
      {
        if (std::isfinite(obs(j))) { seen(n_seen++) = j; }
      }
      seen.resize(n_seen);
      if (n_seen < q)
      {
        obs = obs.elem(seen);
        design = design.rows(seen);
        meas_var = meas_var.submat(seen, seen);
      }

      r_cov.slice(t) = arma::symmatu(r_cov.slice(t));
      const arma::mat rz = r_cov.slice(t) * design.t(); // n_s x n_seen
      const arma::mat f = design * rz + meas_var;
      const arma::vec e = obs - design * s_pred.col(t);

      if (pass == passes - 1)
      {
        loglik(t) = log_density(e, f);
      }

      const arma::mat gain = right_solve(rz, f); // n_s x q

      s_upd.col(t) = s_pred.col(t) + gain * e;
      s_cov.slice(t) = arma::symmatu(r_cov.slice(t) - gain * rz.t());
    }

    // Backward pass. Unlike the two coefficient blocks the transition is not the
    // identity, so this is the full Rauch-Tung-Striebel recursion.
    //
    // A deviation from the code published with the paper, which never fills the
    // lag blocks of the state and so predicts from B_1 alone. That is invisible
    // at p = 1 and wrong above it, so the blocks are carried here.
    arma::mat s_new = s_upd;
    for (arma::uword t = tt - 1; t-- > 0;)
    {
      const arma::mat gain =
        right_solve(s_cov.slice(t) * comp.slice(t).t(), r_cov.slice(t + 1));
      s_new.col(t) += gain * (s_new.col(t + 1) - s_pred.col(t + 1));
    }

    factors = s_new.rows(n_y, n_z - 1).t();
  }

  Rcpp::List posterior =
    Rcpp::List::create(Rcpp::Named("factors") = factors,
                       Rcpp::Named("lambda") = lambda,
                       Rcpp::Named("beta") = beta,
                       Rcpp::Named("u_sigma") = u_sigma,
                       Rcpp::Named("v_sigma") = v_sigma,
                       Rcpp::Named("loglik") = loglik);

  object["posterior"] = posterior;
  return object;
}
