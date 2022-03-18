#' @export
monanova <- function(formula, data, weights, method = "kruskal",
                     verbose = TRUE, max_iter = 1000,
                     rel_tol = 1e-5, degree = 2,
                     difference = 2, knots = 2:10, k = -1,
                     lambda_range = c(0, 1), ...) {

  # data must contain the response (y) and the predictors in the formula
  # k is the number of knots. If k = -1, then it is automatically selected.

  y <- standardize(data$y)
  X <- model.matrix(formula, data = data)
  p <- ncol(X)
  if(is.null(weights)) weights <- rep(1, length(y))

  if(method == "kruskal") {

    fit <- kruskal_monanova(y, X, weights, verbose, max_iter, rel_tol)
    return(fit)

  } else if(method == "smooth1") { # only one step

    fit <- kruskal_monanova(y, X, weights, verbose, max_iter, rel_tol)
    eta <- fit$eta

    boundaries <- range(y) + c(-0.5, 0.5)
    g <- GCV_cpp(eta, y, knots, degree, difference,
                  boundaries, lambda_range)
    k <- g$k
    spline <- g$spline
    D <- g$D
    lambda <- g$lambda
    repeat{
      monotonic_fit <- monotonic(spline, eta, weights, D, lambda)
      if(!anyNA(monotonic_fit$fitted_values)) break
    }
    transformed_y <- standardize(monotonic_fit$fitted_values)
    # if(anyNA(transformed_y)) {
    #   message("NAs encountered")
    #   return(list(spline, eta, weights, D, lambda))
    # }
    linear_model <- lm(transformed_y ~ 0 + X, weights = weights)
    pars <- matrix(linear_model$coefficients)
    sigma2_hat <- summary(linear_model)$sigma^2

    result <- list()
    result$stress <- stress(transformed_y, eta, weights)
    result$sigma2_hat <- sigma2_hat
    result$pars <- pars
    result$eta <- eta
    result$transformed_y <- transformed_y

    result$kruskal <- fit
    result$AIC <- AIC(linear_model)
    result$BIC <- BIC(linear_model)
    result$GCV <- g
    # result$smooth <- monotonic_fit
    result$input <- list(pars = pars, formula = formula, data = data,
                         weights = weights, method = method,
                         verbose = verbose,  max_iter = max_iter,
                         rel_tol = rel_tol, degree = degree,
                         difference = difference, knots = knots, k = k,
                         lambda_range = lambda_range, scam_arguments = ...)

    return(result)

    } else if(method == "smooth2") { # Many steps

      boundaries <- range(y) + c(-0.5, 0.5)
      x_new <- 0

      for(i in 1:max_iter) {

        x_old <- x_new
        opt <- stress_optim(y, X, weights)
        eta <- X %*% opt$pars
        g <- GCV_cpp(eta, y, knots, degree, difference,
                     boundaries, lambda_range)
        transformed_y <- g$transformed_y
        x_new <- stress(transformed_y, eta, weights)

        if(verbose) cat("Iteration ",  i, ": ", x_new, "\n")

        if(abs(x_old - x_new) < rel_tol) break

      }

      k <- g$k
      spline <- g$spline
      D <- g$D
      lambda <- g$lambda
      repeat{
        monotonic_fit <- monotonic(spline, eta, weights, D, lambda)
        if(!anyNA(monotonic_fit$fitted_values)) break
      }
      transformed_y <- standardize(monotonic_fit$fitted_values)

      linear_model <- lm(transformed_y ~ 0 + X, weights = weights)
      pars <- matrix(linear_model$coefficients)
      sigma2_hat <- summary(linear_model)$sigma^2

      result <- list()
      result$AIC <- AIC(linear_model)
      result$BIC <- BIC(linear_model)
      result$stress <- stress(transformed_y, eta, weights)
      result$sigma2_hat <- sigma2_hat
      result$pars <- pars
      result$eta <- eta
      result$transformed_y <- transformed_y

      # result$smooth <- monotonic_fit
      result$input <- list(pars = pars, formula = formula, data = data,
                           weights = weights, method = method,
                           verbose = verbose,  max_iter = max_iter,
                           rel_tol = rel_tol, degree = degree,
                           difference = difference, knots = knots, k = k,
                           lambda_range = lambda_range, scam_arguments = ...)

      return(result)

    } else if(method == "scam1") {

    fit <- kruskal_monanova(y, X, weights, verbose, max_iter, rel_tol)
    eta <- fit$eta
    monotonic_fit <- scam(eta ~ s(y, k = k, bs = "mpi"), weights = weights, ...)

    transformed_y <- standardize(monotonic_fit$fitted.values)
    # lm.wfit(X, y, weight, offset = NULL, method = "qr", tol = 1e-7)
    linear_model <- lm(transformed_y ~ 0 + X, weights = weights)
    pars <- matrix(linear_model$coefficients)
    sigma2_hat <- summary(linear_model)$sigma^2

    result <- list()
    result$stress <- stress(transformed_y, eta, weights)
    result$sigma2_hat <- sigma2_hat
    result$pars <- pars
    result$eta <- eta
    result$transformed_y <- transformed_y

    result$AIC <- AIC(monotonic_fit)
    result$BIC <- BIC(monotonic_fit)
    result$kruskal <- fit
    result$smooth <- monotonic_fit
    result$input <- list(pars = pars, formula = formula, data = data,
                         weights = weights, method = method,
                         verbose = verbose,  max_iter = max_iter,
                         rel_tol = rel_tol, k = k, scam_arguments = ...)

    return(result)

  } else if(method == "scam2") {

    x <- 0

    for(i in 1:max_iter) {

      x_old <- x
      opt <- stress_optim(y, X, weights)
      eta <- X %*% opt$par
      monotonic_fit <- scam(eta ~ s(y, k = k, bs = "mpi"), ...)
      transformed_y <- monotonic_fit$fitted.values
      x <- stress(transformed_y, eta, weights)

      if(verbose) cat("Iteration ",  i, ": ", x, "\n")

      if(abs(x_old - x) < rel_tol) break

    }

    transformed_y <- standardize(monotonic_fit$fitted.values)
    # lm.wfit(X, y, weight, offset = NULL, method = "qr", tol = 1e-7)
    linear_model <- lm(transformed_y ~ 0 + X, weights = weights)
    pars <- matrix(linear_model$coefficients)
    sigma2_hat <- summary(linear_model)$sigma^2

    result <- list()
    result$stress <- stress(transformed_y, eta, weights)
    result$sigma2_hat <- sigma2_hat
    result$pars <- pars
    result$eta <- eta
    result$transformed_y <- transformed_y

    result$AIC <- AIC(monotonic_fit)
    result$BIC <- BIC(monotonic_fit)
    result$smooth <- monotonic_fit
    result$input <- list(pars = pars, formula = formula, data = data,
                         weights = weights, method = method,
                         verbose = verbose,  max_iter = max_iter,
                         rel_tol = rel_tol, k = k, scam_arguments = ...)

    return(result)

  } else {
    stop("Available methods: c('kruskal', 'smooth', 'scam')")
  }

}

