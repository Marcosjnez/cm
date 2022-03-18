// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <ensmallen.hpp>

typedef std::vector<double> stdvec;

double stress(arma::mat beta, arma::vec y, arma::mat X, arma::vec w) {

  arma::mat W = arma::diagmat(w);
  arma::vec eta = X * beta;
  double mu_mean = arma::as_scalar(eta.t() * w) / sum(w);
  arma::vec mu_centered = eta - mu_mean;
  arma::vec residuals = y - eta;
  double u = arma::as_scalar( residuals.t() * W * residuals );
  double v = arma::as_scalar( mu_centered.t() * W * mu_centered );
  double stress = u/v;

  return stress;

}

// [[Rcpp::export]]
arma::vec standardize(arma::vec x){

  x = (x - arma::mean(x)) / arma::stddev(x);

  return x;
}

std::vector<int> sort_indexes(stdvec v) {

  std::vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  stable_sort(idx.begin(), idx.end(), [&v](int i1, int i2) {
    return v[i1] < v[i2];
  });

  return idx;
}

arma::vec int_pava(arma::vec a, arma::uvec c, arma::vec b) {

  stdvec y = arma::conv_to< stdvec >::from(a);
  stdvec w = arma::conv_to< stdvec >::from(b);
  stdvec order = arma::conv_to< stdvec >::from(c);
  int n = y.size();
  int max_idx;
  //int k = 1;
  stdvec result(n);
  stdvec yhat(n);
  stdvec what(n);
  stdvec diffs(n-1);
  std::vector<int> idx(n);
  stdvec wyhat(n);
  std::vector<int> index(n);
  stdvec sumwyhat(n);
  stdvec new_w(n);
  stdvec new_yhat(n);
  // w and y according to the monotonic ordering of x:
  for(int i = 0; i < n; i++) {
    yhat[i] = y[order[i]];
    //Rcout << "yhat is " << yhat[i] << "\n";
    what[i] = w[order[i]];
    //Rcout << "what is " << what[i] << "\n";
    if(i > 0) diffs[i-1] = yhat[i] - yhat[i-1];
    //Rcout << "diffs is " << diffs[i-1] << "\n";
  }
  while(*min_element(diffs.begin(), diffs.end()) < 0) {
    //Rcout << "iteration " << k << "\n";
    //k++;
    for(int i = 0; i < n; i++) {
      wyhat[i] = what[i] * yhat[i];
      if(i > 0) idx[i] = idx[i-1] + (diffs[i-1] > 0);
      //Rcout << "idx is " << idx[i] << "\n";
    }
    max_idx = *max_element(idx.begin(), idx.end());
    sumwyhat.resize(max_idx);
    new_w.resize(max_idx);
    new_yhat.resize(max_idx);
    for(int j = 0; j <= max_idx; ++j) {
      index[j] = j;
      sumwyhat[j] = 0;
      new_w[j] = 0;
      new_yhat[j] = 0;
      //Rcout << "index is " << index[j] << "\n";
    }
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j <= max_idx; ++j) {
        if(idx[i] == index[j]) {
          sumwyhat[j] = (sumwyhat[j] + wyhat[i]);
          new_w[j] = (new_w[j] + what[i]);
        }
      }
    }
    for(int j = 0; j <= max_idx; j++) {
      new_yhat[j] = sumwyhat[j] / new_w[j];
      //Rcout << "new_yhat is " << new_yhat[j] << "\n";
    }
    for(int i = 0; i < n; ++i) {
      yhat[i] = new_yhat[idx[i]];
      //Rcout << "yhat is " << yhat[i] << "\n";
      if(i > 0) diffs[i-1] = yhat[i] - yhat[i-1];
      //Rcout << "diffs is " << diffs[i-1] << "\n";
    }
  }
  for(int i = 0; i < n; ++i) {
    result[order[i]] = yhat[i];
    //Rcout << "result is " << result[i] << "\n";
  }
  arma::vec result2 = arma::conv_to< arma::vec >::from(result);

  return result2;

}

// [[Rcpp::export]]
arma::vec pava(arma::vec a, arma::vec c, arma::vec b) {

  stdvec y = arma::conv_to< stdvec >::from(a);
  stdvec w = arma::conv_to< stdvec >::from(b);
  arma::uvec v = arma::sort_index(c);
  stdvec order = arma::conv_to< stdvec >::from(v);
  int n = y.size();
  int max_idx;
  //int k = 1;
  stdvec result(n);
  stdvec yhat(n);
  stdvec what(n);
  stdvec diffs(n-1);
  std::vector<int> idx(n);
  stdvec wyhat(n);
  std::vector<int> index(n);
  stdvec sumwyhat(n);
  stdvec new_w(n);
  stdvec new_yhat(n);
  // w and y according to the monotonic ordering of x:
  for(int i = 0; i < n; i++) {
    yhat[i] = y[order[i]];
    //Rcout << "yhat is " << yhat[i] << "\n";
    what[i] = w[order[i]];
    //Rcout << "what is " << what[i] << "\n";
    if(i > 0) diffs[i-1] = yhat[i] - yhat[i-1];
    //Rcout << "diffs is " << diffs[i-1] << "\n";
  }
  while(*min_element(diffs.begin(), diffs.end()) < 0) {
    //Rcout << "iteration " << k << "\n";
    //k++;
    for(int i = 0; i < n; i++) {
      wyhat[i] = what[i] * yhat[i];
      if(i > 0) idx[i] = idx[i-1] + (diffs[i-1] > 0);
      //Rcout << "idx is " << idx[i] << "\n";
    }
    max_idx = *max_element(idx.begin(), idx.end());
    sumwyhat.resize(max_idx);
    new_w.resize(max_idx);
    new_yhat.resize(max_idx);
    for(int j = 0; j <= max_idx; ++j) {
      index[j] = j;
      sumwyhat[j] = 0;
      new_w[j] = 0;
      new_yhat[j] = 0;
      //Rcout << "index is " << index[j] << "\n";
    }
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j <= max_idx; ++j) {
        if(idx[i] == index[j]) {
          sumwyhat[j] = (sumwyhat[j] + wyhat[i]);
          new_w[j] = (new_w[j] + what[i]);
        }
      }
    }
    for(int j = 0; j <= max_idx; j++) {
      new_yhat[j] = sumwyhat[j] / new_w[j];
      //Rcout << "new_yhat is " << new_yhat[j] << "\n";
    }
    for(int i = 0; i < n; ++i) {
      yhat[i] = new_yhat[idx[i]];
      //Rcout << "yhat is " << yhat[i] << "\n";
      if(i > 0) diffs[i-1] = yhat[i] - yhat[i-1];
      //Rcout << "diffs is " << diffs[i-1] << "\n";
    }
  }
  for(int i = 0; i < n; ++i) {
    result[order[i]] = yhat[i];
    //Rcout << "result is " << result[i] << "\n";
  }
  arma::vec result2 = arma::conv_to< arma::vec >::from(result);

  return result2;

}

class argmin_stress
{
public:

  argmin_stress(arma::vec y, arma::mat X, arma::vec w) :
  y(y), X(X), w(w) { }

  double EvaluateWithGradient(const arma::mat& par, arma::mat& gradient)
  {
    arma::mat W = arma::diagmat(w);
    arma::vec eta = X * par;
    double mu_mean = arma::as_scalar(eta.t() * w) / sum(w);
    arma::vec mu_centered = eta - mu_mean;
    arma::vec residuals = y - eta;
    double u = arma::as_scalar( residuals.t() * W * residuals );
    double v = arma::as_scalar( mu_centered.t() * W * mu_centered );
    double f = u/v;

    gradient = -2 * X.t() * W * (residuals/v + mu_centered * u / (v*v));
    return f;

  }

private:
  arma::vec y;
  arma::mat X;
  arma::vec w;
};

std::tuple<double, arma::mat> stress_opt(arma::mat par, arma::vec y, arma::mat X, arma::vec w,
                                         int numBasis = 10, int maxIterations = 1000,
                                         double armijoConstant = 1e-4, double wolfe = 1,
                                         double minGradientNorm = 1e-6, double factr = 1e-8,
                                         int maxLineSearchTrials = 50, double minStep = 1e-15,
                                         double maxStep = 1e20) {

  argmin_stress rf(y, X, w);
  ens::L_BFGS lbfgs = ens::L_BFGS(numBasis, maxIterations, armijoConstant,
                                  wolfe, minGradientNorm, factr,
                                  maxLineSearchTrials, minStep, maxStep);

  lbfgs.Optimize(rf, par);
  double f = stress(par, y, X, w);

  std::tuple<double, arma::mat> result = std::make_tuple(f, par);

  return result;

}

// [[Rcpp::export]]
Rcpp::List stress_optim(arma::vec y, arma::mat X, arma::vec weights,
                        int numBasis = 10, int maxIterations = 1000,
                        double armijoConstant = 1e-4, double wolfe = 1,
                        double minGradientNorm = 1e-6, double factr = 1e-8,
                        int maxLineSearchTrials = 50, double minStep = 1e-15,
                        double maxStep = 1e20) {

  Rcpp::List result;

  arma::mat pars(X.n_cols, 1, arma::fill::randu);
  argmin_stress rf(y, X, weights);
  ens::L_BFGS lbfgs = ens::L_BFGS(numBasis, maxIterations, armijoConstant,
                                  wolfe, minGradientNorm, factr,
                                  maxLineSearchTrials, minStep, maxStep);

  lbfgs.Optimize(rf, pars);
  double f = stress(pars, y, X, weights);

  result["stress"] = f;
  result["pars"] = pars;

  return result;

}

// [[Rcpp::export]]
Rcpp::List kruskal_monanova(arma::vec y, arma::mat X, arma::vec w,
                            bool verbose = false, int max_iter = 100,
                            double rel_tol = 1e-5) {

  Rcpp::Timer timer;
  Rcpp::List input, info, result;

  arma::mat pars(X.n_cols, 1, arma::fill::randu);
  // Parameters to pass to the Ensmallen function Ens::L_BFGS:
  int numBasis = 10;
  int maxIterations = 1000;
  double armijoConstant = 1e-4;
  double wolfe = 1;
  double minGradientNorm = 1e-6;
  double factr = 1e-8;
  int maxLineSearchTrials = 50;
  double minStep = 1e-15;
  double maxStep = 1e20;

  int n = y.size();
  int p = X.n_cols;
  double std_y = arma::stddev(y);
  double mean_y = arma::mean(y);

  arma::uvec ordering = arma::sort_index(y);
  double stress = 0;
  arma::vec stress_vector(max_iter);
  double old_stress;
  arma::vec eta;

  int i = 0;

  std::tuple<double, arma::mat> fit;

  do{

    old_stress = stress;

    fit = stress_opt(pars, y, X, w,
                     numBasis, maxIterations,
                     armijoConstant, wolfe,
                     minGradientNorm, factr,
                     maxLineSearchTrials, minStep,
                     maxStep);

    stress = std::get<0>(fit);
    pars = std::get<1>(fit);
    stress_vector[i] = stress;

    eta = X * pars;
    y = int_pava(eta, ordering, w);
    // arma::vec ordering2 = arma::conv_to<arma::vec>::from(ordering);
    // y = pava(eta, ordering2, w);

    ++i;

    if(verbose) {
      Rcpp::Rcout << "Iteration " << i << ": " << stress << "\n";
    }

  } while (std::abs(stress - old_stress) > rel_tol && i < max_iter);

  bool convergence = std::abs(stress - old_stress) < rel_tol && i <= max_iter;

  if(!convergence) {
    Rcpp::warning("Failed convergence in monanova \n Consider either increasing the number of iterations or decreasing the relative tolerance");
  }

  info["stress"] = stress_vector(arma::span(0, i-1));

  input["y"] = y;
  input["X"] = X;
  input["w"] = w;
  input["max_iter"] = max_iter;
  input["rel_tol"] = rel_tol;

  result["stress"] = stress;
  result["sigma2_hat"] = stress*n / (n-p);
  result["pars"] = pars;
  result["eta"] = standardize(eta);
  result["transformed_y"] = standardize(y);
  if(convergence) {
    result["convergence"] = "yes";
  } else {
    result["convergence"] = "no";
  }
  result["iterations"] = i;
  result["info"] = info;
  result["input"] = input;

  timer.step("elapsed");
  result["elapsed"] = timer;

  result.attr("class") = "monanova";

  return result;

}

class penalized_monotonic_regression
{
public:

  penalized_monotonic_regression(arma::mat X, arma::vec y, arma::vec w,
                                 arma::mat S, double lambda) :
  X(X), y(y), w(w), S(S), lambda(lambda) { }

  double EvaluateWithGradient(const arma::mat& beta, arma::mat& gradient)
  {
    int p = beta.n_rows;
    arma::mat W = arma::diagmat(w);
    arma::mat trans_beta = beta;
    for(int i=1; i < p; i++) trans_beta[i] = exp(beta[i]);
    arma::vec eta = X * trans_beta;
    arma::vec residuals = y - eta;
    double objective_value = arma::as_scalar( residuals.t() * W * residuals ) +
      lambda * arma::as_scalar( beta.t() * S * beta);

    gradient = -2 * X.t() * W * residuals;
    for(int i=1; i < p; i++) gradient[i] *= trans_beta[i];
    gradient += 2 * lambda * S * beta;

    return objective_value;

  }

private:
  arma::mat X;
  arma::vec y;
  arma::vec w;
  arma::mat S;
  double lambda;
};

// [[Rcpp::export]]
Rcpp::List monotonic(arma::mat X, arma::vec y, arma::vec weights,
                     arma::mat D, double lambda,
                     int numBasis = 10, int maxIterations = 1000,
                     double armijoConstant = 1e-4, double wolfe = 2,
                     double minGradientNorm = 1e-6, double factr = 1e-8,
                     int maxLineSearchTrials = 50, double minStep = 1e-15,
                     double maxStep = 1e20) {

  arma::mat beta(X.n_cols, 1, arma::fill::randu);
  arma::mat S = D.t() * D;
  penalized_monotonic_regression rf(X, y, weights, S, lambda);

  ens::L_BFGS lbfgs = ens::L_BFGS(numBasis, maxIterations, armijoConstant,
                                  wolfe, minGradientNorm, factr,
                                  maxLineSearchTrials, minStep, maxStep);
  lbfgs.Optimize(rf, beta);

  arma::mat trans_beta = beta;
  for(int i=1; i < beta.size(); i++) trans_beta[i] = exp(beta[i]);
  // arma::vec fitted_values = X * trans_beta;
  arma::mat fitted_values = X * trans_beta;

  Rcpp::List result;
  result["pars"] = beta;
  result["fitted_values"] = fitted_values;

  Rcpp::List input = Rcpp::List::create(Rcpp::_("X") = X, Rcpp::_["y"] = y,
                                        Rcpp::_["weights"] = weights,
                                        Rcpp::_["D"] = D, Rcpp::_["lambda"] = lambda);
  result["input"] = input;

  return result;
}

arma::mat center_mat(arma::mat X, int n, int p) {

  arma::rowvec means(p);
  double sums;

  for(int j = 0; j < p; j++) {
    sums = 0;
    for(int i = 0; i < n; i++) {
      sums += X(i, j);
    }
    means(j) = sums/n;
  }

  X.each_row() -= means;

  return X;
}

// [[Rcpp::export]]
arma::mat b_spline(arma::vec x, arma::vec knots, int degree, arma::vec boundaries,
                   bool center = false, bool intercept = false) {

  int n = x.size();
  unsigned int p = knots.size() + degree;
  double z0, z1, output;

  arma::vec lower_boundary(degree);
  arma::vec upper_boundary(degree);
  for(int i = 0; i < degree; i++) {
    lower_boundary(i) = boundaries(0);
  }
  for(int i = 0; i < degree; i++) {
    upper_boundary(i) = boundaries(1);
  }

  knots = arma::join_cols(lower_boundary, knots, upper_boundary);

  arma::mat X(n, p);

  for(int j = 0; j < n; j++) {
    for(int k = 0; k < p; k++) {
      if((x[j] <= knots[k+1]) && (x[j] > knots[k])) {
        X(j, k) = 1;
      } else {
        X(j, k) = 0;
      }
    }
  }

  X.insert_cols(p, 1);

  for(int m = 1; m < (degree+1); m++) {
    for(int j = 0; j < n; j++) {
      for(int k = 0; k < p; k++) {

        if(knots[k+m] == knots[k]) {
          z0 = 0;
        } else {
          z0 = (x[j] - knots[k]) / (knots[k+m] - knots[k]);
        }
        if(knots[k+m+1] == knots[k+1]) {
          z1 = 0;
        } else {
          z1 = (knots[k+m+1] - x[j]) / (knots[k+m+1] - knots[k+1]);
        }

        X(j, k) = z0 * X(j, k) + z1 * X(j, k+1);
      }
    }
  }

  arma::uvec indices = {p};
  X.shed_cols(indices);

  if(center) {
    X = center_mat(X, n, p);
  }

  if(intercept) {
    X.insert_cols(0, arma::ones(n));
  }

  return X;
}

// [[Rcpp::export]]
arma::mat tri_mat(int dim) {

  arma::mat matrix(dim, dim, arma::fill::zeros);
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      if(i >= j) matrix(i, j) = 1;
    }
  }

  return matrix;
}

// [[Rcpp::export]]
arma::mat diff_mat(int n, int p, int difference, bool intercept = true,
                   bool multiply = true) { // false

  arma::mat X = arma::eye(n, p);
  X = diff(X, difference);

  if(multiply) {

    double prod = 1 + 0.5*(p-3);
    for(int i=0; i < X.n_rows; ++i) {
      for(int j=0; j < X.n_cols; ++j) {
        X(i, j) *= prod;
      }
    }

  }

  if(intercept) X.insert_cols(0, 1);

  return X;
}

arma::vec seq_cpp(double from, double to, int length_out) {

  double interval = (to - from) / (length_out-1);
  arma::vec sequence(length_out);
  for(int i=0; i < length_out; i++) sequence(i) = from + interval*i;

  return sequence;

}

double GCV_score(arma::vec y, arma::mat spline,
                 arma::mat D, double lambda) {

  int n = spline.n_rows;

  arma::mat lambda_D = sqrt(lambda) * D;
  arma::mat X = arma::join_cols(spline, lambda_D);

  arma::mat XX = X.t() * X;
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::svd(U, s, V, XX, "dc");
  arma::mat inv_XX = U * arma::diagmat(1/s) * V.t();
  arma::mat H = X * inv_XX * X.t();
  arma::vec h = arma::diagvec(H);
  arma::vec residuals = y - H*y;
  arma::vec residuals2 = residuals % residuals;

  double numerator = n * arma::accu(residuals2(arma::span(0, n-1)));
  double denominator = n - 1.5*arma::accu(h(arma::span(0, n-1))); // 1.5
  denominator *= denominator;
  double score = numerator / denominator;

  return score;

}

Rcpp::List GCV_cpp_old(arma::vec y, arma::vec x, arma::vec lambdas,
                       arma::vec k, int degree,
                       int difference, arma::vec boundaries) {

  //   int cores = omp_get_max_threads();
  //   omp_set_num_threads(cores);
  // #pragma omp parallel for

  int n = y.size();
  int N = lambdas.size();
  arma::mat spline;
  int K = k.size();

  arma::mat Vg(N, K);

  x = standardize(x);
  double min_x = boundaries(0);
  double max_x = boundaries(1);

  for(int a=0; a < K; a++) {

    int p = k[a] + degree + 1;
    arma::mat D = diff_mat(p-1, p-1, difference);
    arma::mat S = D.t() * D;

    arma::vec knots = seq_cpp(min_x, max_x, k[a]+2);
    arma::vec inner_knots = knots(arma::span(1, k[a]));
    spline = b_spline(x, inner_knots, degree, boundaries, false, false) * tri_mat(p-1);
    spline.insert_cols(0, arma::ones(n));

    unsigned int size_aug = n + D.n_rows;
    arma::mat y_aug_temp = arma::conv_to<arma::mat>::from(y);;
    y_aug_temp.resize(size_aug, 1);
    arma::vec y_aug = arma::conv_to<arma::vec>::from(y_aug_temp);

    for(int i = 0; i < N; ++i) {

      Vg(i, a) = GCV_score(y_aug, spline, D, lambdas[i]);

    }

  }

  arma::vec Vg_mins(K);
  for(int i=0; i < K; i++) Vg_mins(i) = min(Vg.col(i));
  arma::uword min_index = index_min(Vg_mins);

  arma::vec selected_Vg = Vg.col(min_index);
  int selected_k = k[min_index];
  double selected_lambda = lambdas[index_min(selected_Vg)];

  Rcpp::List ret;

  ret["k"] = selected_k;
  ret["lambda"] = selected_lambda;
  ret["min_GCV"] = min(Vg_mins);
  ret["lambdas_GCV"] = join_rows(lambdas, selected_Vg);
  ret["Vg"] = Vg;

  ret.attr("class") = "GCV";

  return ret;

}

// [[Rcpp::export]]
Rcpp::List GCV_cpp(arma::vec y, arma::vec x, arma::vec k,
                   int degree, int difference,
                   arma::vec boundaries, arma::vec lambda_range) {

  //   int cores = omp_get_max_threads();
  //   omp_set_num_threads(cores);
  // #pragma omp parallel for

  int n = y.size();
  int K = k.size();

  double numerator;
  double denominator;

  x = standardize(x);
  double min_x = boundaries(0);
  double max_x = boundaries(1);
  arma::mat Vq(K, 4);

  for(int j=0; j < K; j++) {

    double lambda, score;
    int iterations = 0;

    int p = k[j] + degree + 1;
    arma::mat D = diff_mat(p-1, p-1, difference);

    arma::vec knots = seq_cpp(min_x, max_x, k[j]+2);
    arma::vec inner_knots = knots(arma::span(1, k[j]));
    arma::mat spline = b_spline(x, inner_knots, degree, boundaries, false, false) * tri_mat(p-1);
    spline.insert_cols(0, arma::ones(n));

    int size_aug = n + D.n_rows;
    arma::vec y_aug(size_aug);
    y_aug(arma::span(0, n-1)) = y;

    arma::vec lambdas = lambda_range;
    arma::vec scores(2);
    scores[0] = GCV_score(y_aug, spline, D, lambdas[0]);
    scores[1] = GCV_score(y_aug, spline, D, lambdas[1]);

    // https://en.wikipedia.org/wiki/Golden-section_search

    double gr = arma::datum::gratio;
    double a = lambdas[0];
    double b = lambdas[1];
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;

    while(std::abs(b - a) > 1e-5) {

      ++iterations;

      double fc = GCV_score(y_aug, spline, D, c);
      double fd = GCV_score(y_aug, spline, D, d);

      if(fc < fd) {
        b = d;
      } else{
        a = c;
      }
      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
    }

    lambda = (b + a) / 2;
    score = GCV_score(y_aug, spline, D, lambda);

    Vq(j, 0) = k[j];
    Vq(j, 1) = lambda;
    Vq(j, 2) = score;
    Vq(j, 3) = iterations;

  }

  arma::uword min_index = index_min(Vq.col(2));

  int selected_k = Vq(min_index, 0);
  double selected_lambda = Vq(min_index, 1);
  double min_score = Vq(min_index, 2);

  int p = selected_k + degree + 1;
  arma::mat D = diff_mat(p-1, p-1, difference);

  arma::vec knots = seq_cpp(min_x, max_x, selected_k+2);
  arma::vec inner_knots = knots(arma::span(1, selected_k));
  arma::mat spline = b_spline(x, inner_knots, degree, boundaries, false, false) * tri_mat(p-1);
  spline.insert_cols(0, arma::ones(n));

  arma::mat lambda_D = sqrt(selected_lambda) * D;
  arma::mat X = arma::join_cols(spline, lambda_D);

  int size_aug = n + D.n_rows;
  arma::vec y_aug(size_aug);
  y_aug(arma::span(0, n-1)) = y;

  arma::mat XX = X.t() * X;
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::svd(U, s, V, XX, "dc");
  arma::mat inv_XX = U * arma::diagmat(1/s) * V.t();
  arma::mat H = X * inv_XX * X.t();
  arma::vec h = arma::diagvec(H);
  arma::vec temp = H*y_aug;
  arma::vec transformed_y = temp(arma::span(0, n-1));

  Rcpp::List result;

  result["transformed_y"] = standardize(transformed_y);
  result["k"] = selected_k;
  result["lambda"] = selected_lambda;
  result["GCV"] = min_score;
  result["GCV_matrix"] = Vq;
  result["spline"] = spline;
  result["D"] = D;

  result.attr("class") = "GCV";

  return result;

}
