# Generic functions:
#' @export
stress <- function(y, eta, w) {

  eta_mean <- sum(eta * w) / sum(w)
  eta_centered <- eta - eta_mean
  residuals <- y - eta
  U <- sum(t(residuals) %*% diag(w) %*% residuals)
  V <- sum(t(eta_centered) %*% diag(w) %*% eta_centered)
  stress <- U/V

  return(stress)

}
#' @export
plot.monanova <- function(fit, ...) {

  x <- fit$y
  y <- fit$transformed.response

  if(fit$type == 2) {
    plot(x/fit$w, y, main = "Best fitting transformation",
         ylab = "Transformed response", xlab = "Response", las = 1, ...)
  } else {
    plot(x, y, main = "Best fitting transformation",
         ylab = "Transformed response", xlab = "Response", las = 1, ...)
  }

}
#' @export
plot.monotonic <- function(fit, ...) {

  if(fit$type == 2) {

    x <- sort(fit$eta)
    y <- (fit$y/fit$w)[order(fit$eta)]

    plot(x, y, main = "Predictions", xlab = "Linear combination of predictors in X",
         ylab = "f(X)", las = 1, type = "n")

    boot <- t(fit$boot_fitted_values)

    intervals <- apply(boot, 2, FUN = quantile, probs = c(.025, .975), na.rm = TRUE)
    intervals <- intervals[, order(fit$eta)]

    polygon(c(rev(x), x), c(rev(intervals[1, ]), intervals[2, ]),
            col = "#fde0e0", border = NA)

    lines(x, intervals[1, ], lty = 'dashed', col = 'black')
    lines(x, intervals[2, ], lty = 'dashed', col = 'black')

    lines(x, fit$fitted.values[order(fit$eta)], lwd = 3, ...)
    points(x, y)

  } else {

    x <- sort(fit$eta)
    y <- fit$y[order(fit$eta)]

    plot(x, y, main = "Predictions", xlab = "Linear combination of predictors in X",
         ylab = "f(X)", las = 1, type = "n")

    if( !is.null(fit$boot_fitted_values) ) {

      boot <- t(fit$boot_fitted_values)

      intervals <- apply(boot, 2, FUN = quantile, probs = c(.025, .975), na.rm = TRUE)
      intervals <- intervals[, order(fit$eta)]

      polygon(c(rev(x), x), c(rev(intervals[1, ]), intervals[2, ]),
              col = "#fde0e0", border = NA)

      lines(x, intervals[1, ], lty = 'dashed', col = 'black')
      lines(x, intervals[2, ], lty = 'dashed', col = 'black')

    }

    lines(sort(x), fit$fitted.values[order(fit$eta)], lwd = 3, ...)
    points(x, y)
  }
}
#' @export
plot.GCV <- function(scores, ...) {

  plot(scores$lambdas_GCV[, 1], scores$lambdas_GCV[, 2], xlab = bquote(log(lambda)), ylab = "", ...)

}

# Functions to generate different composition rules:
#' @export
additive <- function(I, J, K) {

  A <- runif(I)
  P <- runif(J)
  U <- runif(K)

  eta <- vector(length = I*J*K)

  index <- 0
  for(i in 1:I) {
    for(j in 1:J) {
      for(k in 1:K) {
        index <- index + 1
        eta[index] <- A[i] + P[j] + U[k]
      }
    }
  }

  eta <- (eta - mean(eta))/sd(eta)

  return(eta)
}
#' @export
distributive <- function(I, J, K) {

  A <- runif(I)
  P <- runif(J)
  U <- runif(K)

  eta <- vector(length = I*J*K)

  index <- 0
  for(i in 1:I) {
    for(j in 1:J) {
      for(k in 1:K) {
        index <- index + 1
        eta[index] <- (A[i] + P[j]) * U[k]
      }
    }
  }

  eta <- (eta - mean(eta))/sd(eta)

  return(eta)
}
#' @export
dual_distributive <- function(I, J, K) {

  A <- runif(I)
  P <- runif(J)
  U <- runif(K)

  eta <- vector(length = I*J*K)

  index <- 0
  for(i in 1:I) {
    for(j in 1:J) {
      for(k in 1:K) {
        index <- index + 1
        eta[index] <- A[i]*U[k] + P[j]
      }
    }
  }

  eta <- (eta - mean(eta))/sd(eta)

  return(eta)
}
#' @export
semi_additive <- function(I, J, K) {

  A <- runif(I)
  P <- runif(J)

  eta <- vector(length = I*J*K)

  index <- 0
  for(i in 1:I) {
    for(j in 1:J) {
      for(k in 1:K) {
        index <- index + 1
        eta[index] <- A[i] + P[j]
      }
    }
  }

  eta <- (eta - mean(eta))/sd(eta)

  return(eta)
}
#' @export
additive_joint <- function(I, J, K) {

  A <- runif(I)
  P <- runif(J)
  U <- runif(K)

  eta <- vector(length = I*J*K)

  index <- 0
  for(i in 1:I) {
    for(j in 1:J) {
      for(k in 1:K) {
        index <- index + 1
        eta[index] <- A[i] + P[j] + U[k] + A[i]*P[j]
      }
    }
  }

  eta <- (eta - mean(eta))/sd(eta)

  return(eta)
}

# Function to generate the design matrix:
#' @export
setup <- function(model, I, J, K, times) {

  fa <- factor(rep(1:I, each = J*K))
  fp <- factor(rep(1:J, times = I, each = K))
  fu <- factor(rep(1:K, times = I*J))

  if(model == "additive") {
    eta <- rep(additive(I, J, K), times)
    X <- model.matrix(~ fa + fp + fu)
  } else if(model == "distributive") {
    eta <- rep(distributive(I, J, K), times)
    X <- model.matrix(~ fa:fu + fp:fu)
  } else if(model == "dual_distributive") {
    eta <- rep(dual_distributive(I, J, K), times)
    X <- model.matrix(~ fa * fu + fp)
  } else if(model == "semi_additive") {
    eta <- rep(additive(I, J, K), times)
    X <- model.matrix(~ fa + fp)
  } else if(model == "additive_joint") {
    eta <- rep(additive(I, J, K), times)
    X <- model.matrix(~ fa * fp + fu)
  }

  n <- nrow(X)
  X <- X[rep(1:n, times), ]
  n <- nrow(X)
  constant <- n / (n - ncol(X))

  return(list(n = n, constant = constant, X = X, eta = eta))

}

# Functions to generate data from a composition rule:
#' @export
sim_data <- function(f, set_up, sigma) {

  eta <- set_up$eta
  n <- set_up$n
  error <- sqrt(sigma^2/(1 - sigma^2))
  y <- f(eta + rnorm(n, 0, error))

  return(y)

}

# Functions to run the simulations:
#' @export
iterations <- function(sigma, S, f, set_up, beta, w,
                       smooth, lambdas, ks) {

  X <- set_up$X
  sigma2_hat <- vector(length = S)

  for(i in 1:S) {

    y <- sim_data(f, set_up, sigma)

    fit <- monanova(beta, y = y, w = w, X = X, type = 1,
                    verbose = FALSE, smooth = smooth,
                    degree = 2, difference = 2,
                    GCV = TRUE, lambdas = lambdas, ks = ks,
                    max_iter = 100, abs_tol = 1e-9, rel_tol = 1e-5)

    sigma2_hat[i] <- fit$stress * set_up$constant

  }

  return(list(sigma2_hat = sigma2_hat))

}
#' @export
simulations <- function(S, set_up, fs, smooth, lambdas, ks, sigmas) {

  n <- set_up$n
  w <- rep(1, n)
  beta <- matrix(runif(ncol(set_up$X)))

  results <- list()
  i <- 0
  for(f in fs) {
    i <- i + 1
    results[[i]] <- mclapply(sigmas, FUN = iterations, S, f,
                             set_up, beta, w, smooth, lambdas, ks,
                             mc.cores = 4L)
  }

  return(results)
}

# Monotonic functions:
#' @export
staggered <- function(eta) {
  x <- 3*eta / max(abs(eta))
  y <- x + sin(2*x) / 2
  y <- y + abs(min(y))
  y <- y/max(y)
  return(y)
}
#' @export
plogis2 <- function(eta) {
  x <- 3*eta / max(abs(eta))
  y <- plogis(2*x)
  return(y)
}
#' @export
cubic <- function(eta) {
  x <- 3*eta / max(abs(eta))
  y <- x^3
  y <- y + abs(min(y))
  y <- y/max(y)
  return(y)
}
#' @export
exponential <- function(eta) {
  x <- 3*eta / max(abs(eta))
  y <- exp(x)
  y <- y + abs(min(y))
  y <- y/max(y)
  return(y)
}
#' @export
logarithmic <- function(eta) {
  x <- eta + abs(min(eta))
  x <- 3*x / max(x)
  y <- log(x + 1e-4) # add 1e-4 in case x=0
  y <- y + abs(min(y))
  y <- y/max(y)
  return(y)
}

# Plot the functions:
# curve(staggered(x), xlim = c(-3, 3), ylab = "y", lwd = 3,
#       las = 1, xlab = "Linear combinations of predictions in X")
# curve(cubic(x), xlim = c(-3, 3), col = "red", lwd = 3, add = TRUE)
# curve(plogis2(x), xlim = c(-3, 3), col = "blue", lwd = 3, add = TRUE)
# curve(exponential(x), xlim = c(-3, 3), col = "green", lwd = 3, add = TRUE)
# x <- seq(-3, 3, length.out = 100)
# lines(x, logarithmic(x), xlim = c(-3, 3), col = "brown", lwd = 3)



