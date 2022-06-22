


internal_R_fitglmm_ai_dense <- function (Y, X, q, kins, ng, group.idx, W, tau, fixtau, k, encode = F, sample_encode = NULL, a = NULL)
{
  n <- nrow(X)
  p <- ncol(X)
  q2 <- sum(fixtau == 0)
  diagSigma <- rep(0, n)
  for (i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/W[group.idx[[i]]]
  
  Sigma <- diag(diagSigma)
  for (i in 1:q) Sigma <- Sigma + tau[i + ng] * kins[[i]]
  
  if (encode){
    a_encode <- a%*%sample_encode
    Sigma_i <- MASS::ginv(((a_encode)%*%Sigma%*%t(a_encode)))
    Y <- a_encode%*%Y
    X <- Y-Y+1
    n <- nrow(X)
    diagSigma <- diagSigma[1:n]
    W <- W[1:n]
    kins[[1]] <- (a_encode)%*%kins[[1]]%*%t(a_encode)
    group.idx <- list(c(1:k))
    
  } else {
    Sigma_i <- (MASS::ginv(Sigma))
  }
  rm(Sigma)
  gc()
  Sigma_iX <- crossprod(Sigma_i, X)
  XSigma_iX <- crossprod(X, Sigma_iX)
  cov <- MASS::ginv(XSigma_iX)
  Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
  P <- Sigma_i - tcrossprod(Sigma_iXcov, Sigma_iX)
  alpha <- crossprod(cov, crossprod(Sigma_iX, Y))
  eta <- Y - diagSigma * (crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX,
                                                             t(alpha)))
  rm(Sigma_i)
  gc()
  
  idxtau <- which(fixtau == 0)
  PY <- crossprod(P, Y)
  wPY <- PY/W
  diagP <- diag(P)/W
  AI <- matrix(NA, q2, q2)
  score <- rep(NA, q2)
  for (i in 1:q2) {
    if (idxtau[i] <= ng) {
      score[i] <- sum(wPY[group.idx[[idxtau[i]]]] *
                        PY[group.idx[[idxtau[i]]]] - diagP[group.idx[[idxtau[i]]]])
      for (j in 1:i) {
        AI[i, j] <- crossprod(wPY[group.idx[[idxtau[i]]]],
                              crossprod(P[group.idx[[idxtau[j]]], group.idx[[idxtau[i]]]],
                                        wPY[group.idx[[idxtau[j]]]]))
        if (j != i)
          AI[j, i] <- AI[i, j]
      }
    }
    else {
      PAPY <- crossprod(P, crossprod(kins[[idxtau[i] -
                                             ng]], PY))
      score[i] <- sum((Y * PAPY)) - sum((P * kins[[idxtau[i] -
                                                     ng]]))
      for (j in 1:i) {
        if (idxtau[j] <= ng)
          AI[i, j] <- sum((wPY[group.idx[[idxtau[j]]]] *
                             PAPY[group.idx[[idxtau[j]]]]))
        else AI[i, j] <- sum((PY * crossprod(kins[[idxtau[j] -
                                                     ng]], PAPY)))
        if (j != i)
          AI[j, i] <- AI[i, j]
      }
    }
  }
  Dtau <- (score%*%t(AI)%*%MASS::ginv(AI%*%t(AI)))
  return(list(y = Y, Dtau = Dtau, P = P, cov = cov, alpha = alpha,
              eta = eta, sample_encode = sample_encode))
}


david_glmmkin.ai <- function (fit0, kins, k, tau = rep(0, length(kins) + length(group.idx)), fixtau = rep(0, length(kins) + length(group.idx)), fixrho = NULL, maxiter = 500, tol = 1e-05, verbose = FALSE, encode = FALSE, sample_encode = NULL) {
  
  
  q <- 1
  ng <- 1
  n.pheno <- ncol(fit0$model[[1]])
  n <- nrow(fit0$model[[1]])
  
  mdl <- model.frame(formula = y~1)
  idx <- match(rownames(mdl), rownames(model.frame(formula = y~1)))
  group.id <- rep(1, length(idx))
  
  group.unique <- unique(group.id)
  group.idx <- list()
  for (i in 1:length(group.unique)) group.idx[[i]] <- which(group.id ==
                                                              group.unique[i])
  
  if (!is.null(n.pheno)) {
    n.params <- n.pheno * (n.pheno + 1)/2 * (q +
                                               ng)
    covariance.idx <- matrix(0, n.pheno * (n.pheno -
                                             1)/2 * (q + ng), 3)
    
    for (ii in 1:ng) {
      diagSigma <- rep(0, n)
      diagSigma[group.idx[[ii]]] <- 1
      for (jj in 1:n.pheno) {
        for (kk in jj:n.pheno) {
          if (kk > jj)
            covariance.idx[(ii - 1) * n.pheno *
                             (n.pheno - 1)/2 + (2 * n.pheno - jj) *
                             (jj - 1)/2 + kk - jj, ] <- c((2 *
                                                             n.pheno + 2 - jj) * (jj - 1)/2 + kk -
                                                            jj + 1, (2 * n.pheno + 2 - jj) * (jj -
                                                                                                1)/2 + 1, (2 * n.pheno + 2 - kk) *
                                                            (kk - 1)/2 + 1) + (ii - 1) * n.pheno *
              (n.pheno + 1)/2
        }
      }
    }
    for (ii in 1:q) {
      for (jj in 1:n.pheno) {
        for (kk in jj:n.pheno) {
          if (kk > jj)
            covariance.idx[(ng + ii - 1) * n.pheno *
                             (n.pheno - 1)/2 + (2 * n.pheno - jj) *
                             (jj - 1)/2 + kk - jj, ] <- c((2 *
                                                             n.pheno + 2 - jj) * (jj - 1)/2 + kk -
                                                            jj + 1, (2 * n.pheno + 2 - jj) * (jj -
                                                                                                1)/2 + 1, (2 * n.pheno + 2 - kk) *
                                                            (kk - 1)/2 + 1) + (ng + ii - 1) *
              n.pheno * (n.pheno + 1)/2
        }
      }
    }
  }
  
  
  
  is.Matrix <- any(sapply(kins, function(xx) !is.null(attr(class(xx),
                                                           "package")) && attr(class(xx), "package") == "Matrix"))
  y <- fit0$y
  n <- length(y)
  offset <- fit0$offset
  if (is.null(offset))
    offset <- rep(0, n)
  family <- gaussian()
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu)/mu.eta
  sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0)) * family$variance(mu))
  X <- model.matrix(fit0)
  alpha <- fit0$coef
  if (verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }
  if (family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }
  q <- length(kins)
  ng <- length(group.idx)
  idxtau <- which(fixtau == 0)
  if (!is.null(covariance.idx))
    idxtau2 <- intersect(covariance.idx[, 1], idxtau)
  q2 <- sum(fixtau == 0)
  if (q2 > 0) {
    tau[idxtau] <- rep(var(Y)/(q + ng), q2)
    if (!is.null(covariance.idx))
      tau[idxtau2] <- 0
    diagSigma <- rep(0, n)
    for (i in 1:ng) diagSigma[group.idx[[i]]] <- tau[i]/sqrtW[group.idx[[i]]]^2
    Sigma <- diag(diagSigma)
    for (i in 1:q) {
      tau[i + ng] <- tau[i + ng]/mean(diag(kins[[i]]))
      Sigma <- Sigma + tau[i + ng] * kins[[i]]
    }
    Sigma_i <- (MASS::ginv(Sigma))
    rm(Sigma, diagSigma)
    gc()
    Sigma_iX <- crossprod(Sigma_i, X)
    XSigma_iX <- crossprod(X, Sigma_iX)
    if (is.Matrix)
      XSigma_iX <- forceSymmetric(XSigma_iX)
    cov <- (MASS::ginv(XSigma_iX))
    Sigma_iXcov <- tcrossprod(Sigma_iX, cov)
    PY <- crossprod(Sigma_i, Y) - tcrossprod(Sigma_iX, t(crossprod(Sigma_iXcov,
                                                                   Y)))
    tau0 <- tau
    
    h2 <- 0.5
    for (i in 1:q2) {
      if (idxtau[i] <= ng)
        tau[idxtau[i]] <- max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 *
                                (sum((PY/sqrtW)[group.idx[[idxtau[i]]]]^2) -
                                   sum(((diag(Sigma_i) - rowSums(Sigma_iX *
                                                                   Sigma_iXcov))/sqrtW^2)[group.idx[[idxtau[i]]]]))/n)
      else {
        APY <- crossprod(kins[[idxtau[i] - ng]], PY)
        PAPY <- crossprod(Sigma_i, APY) - tcrossprod(Sigma_iX,
                                                     t(crossprod(Sigma_iXcov, APY)))
        if (!is.null(covariance.idx))
          tau[idxtau[i]] <- if (idxtau[i] %in% idxtau2)
            0
        else max(0, tau0[idxtau[i]] + tau0[idxtau[i]]^2 *
                   (sum(Y * PAPY) - (sum(Sigma_i * kins[[idxtau[i] -
                                                           ng]]) - sum(Sigma_iX * crossprod(kins[[idxtau[i] -
                                                                                                    ng]], Sigma_iXcov))))/n)
        else tau[idxtau[i]] <- max(0, tau0[idxtau[i]] +
                                     tau0[idxtau[i]]^2 * (sum(Y * PAPY) - (sum(Sigma_i *
                                                                                 kins[[idxtau[i] - ng]]) - sum(Sigma_iX *
                                                                                                                 crossprod(kins[[idxtau[i] - ng]], Sigma_iXcov))))/n)
      }
    }
  }
  
  
  
  
  
  
  
  
  
  sample_encode <- if(is.null(sample_encode)){diag(n)}else{sample_encode}
  
  if (encode){
    a <- MASS::ginv(sample_encode%*%t(sample_encode))
  }
  
  
  
  for (i in seq_len(maxiter)) {
    if (verbose)
      cat("\nIteration ", i, ":\n")
    alpha0 <- alpha
    tau0 <- tau
    
    fit <- internal_R_fitglmm_ai_dense(Y = Y, X = X, q = q, kins = kins, ng = ng, group.idx = group.idx, W = sqrtW^2, tau = tau, fixtau = fixtau, k=k, encode = encode, sample_encode = sample_encode, a = a)
    
    if (q2 > 0) {
      Dtau <- as.numeric(fit$Dtau)
      tau[idxtau] <- tau0[idxtau] + Dtau
      if (is.null(covariance.idx)) {
        tau[tau < tol & tau0 < tol] <- 0
        while (any(tau < 0)) {
          Dtau <- Dtau/2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[tau < tol & tau0 < tol] <- 0
        }
        tau[tau < tol] <- 0
      }
      else {
        fixrho.idx0 <- apply(covariance.idx, 1, function(x) abs(tau0[x[1]]) >
                               (1 - 1.01 * tol) * sqrt(tau0[x[2]] * tau0[x[3]]))
        tau[-covariance.idx[, 1]][tau[-covariance.idx[,
                                                      1]] < tol & tau0[-covariance.idx[, 1]] < tol] <- 0
        if (any(fixrho != 0)) {
          idxrho <- which(fixrho != 0)
          tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] *
                                                               sqrt(tau[covariance.idx[idxrho, 2]] * tau[covariance.idx[idxrho,
                                                                                                                        3]]))
        }
        fixrho.idx <- suppressWarnings(apply(covariance.idx,
                                             1, function(x) abs(tau[x[1]]) > (1 - 1.01 *
                                                                                tol) * sqrt(tau[x[2]] * tau[x[3]])))
        tau[covariance.idx[fixrho.idx & fixrho.idx0,
                           1]] <- sign(tau[covariance.idx[fixrho.idx &
                                                            fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx &
                                                                                                          fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx &
                                                                                                                                                  fixrho.idx0, 3]])
        while (any(tau[-covariance.idx[, 1]] < 0) ||
               any(apply(covariance.idx, 1, function(x) abs(tau[x[1]]) >
                         sqrt(tau[x[2]] * tau[x[3]])))) {
          Dtau <- Dtau/2
          tau[idxtau] <- tau0[idxtau] + Dtau
          tau[-covariance.idx[, 1]][tau[-covariance.idx[,
                                                        1]] < tol & tau0[-covariance.idx[, 1]] <
                                      tol] <- 0
          if (any(fixrho != 0)) {
            idxrho <- which(fixrho != 0)
            tau[covariance.idx[idxrho, 1]] <- suppressWarnings(fixrho[idxrho] *
                                                                 sqrt(tau[covariance.idx[idxrho, 2]] *
                                                                        tau[covariance.idx[idxrho, 3]]))
          }
          fixrho.idx <- suppressWarnings(apply(covariance.idx,
                                               1, function(x) abs(tau[x[1]]) > (1 - 1.01 *
                                                                                  tol) * sqrt(tau[x[2]] * tau[x[3]])))
          tau[covariance.idx[fixrho.idx & fixrho.idx0,
                             1]] <- sign(tau[covariance.idx[fixrho.idx &
                                                              fixrho.idx0, 1]]) * sqrt(tau[covariance.idx[fixrho.idx &
                                                                                                            fixrho.idx0, 2]] * tau[covariance.idx[fixrho.idx &
                                                                                                                                                    fixrho.idx0, 3]])
        }
        tau[-covariance.idx[, 1]][tau[-covariance.idx[,
                                                      1]] < tol] <- 0
        tau[covariance.idx[fixrho.idx, 1]] <- sign(tau[covariance.idx[fixrho.idx,
                                                                      1]]) * sqrt(tau[covariance.idx[fixrho.idx,
                                                                                                     2]] * tau[covariance.idx[fixrho.idx, 3]])
      }
    }
    cov <- as.matrix(fit$cov)
    alpha <- as.numeric(fit$alpha)
    eta <- if(encode){as.numeric(t(sample_encode)%*%fit$eta) + offset}else{(fit$eta) + offset}
    if (verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu)/mu.eta
    sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0)[1:dim(sample_encode)[1]])*family$variance(mu))
    if (2 * max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) +
                                     tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) <
        tol)
      break
    if (max(abs(tau)) > tol^(-2)) {
      warning("Large variance estimate observed in the iterations, model not converged...",
              call. = FALSE)
      i <- maxiter
      break
    }
  }
  converged <- ifelse(i < maxiter, TRUE, FALSE)
  if (!is.Matrix)
    fit$Sigma_i <- fit$Sigma_iX <- NULL
  names(alpha) <- names(fit0$coef)
  rownames(cov) <- rownames(vcov(fit0))
  colnames(cov) <- colnames(vcov(fit0))
  return(list(theta = tau, n.pheno = 1, n.groups = ng, coefficients = alpha,
              linear.predictors = eta, fitted.values = mu, Y = Y,
              X = X, P = fit$P, cov = cov, Sigma_i = fit$Sigma_i,
              Sigma_iX = fit$Sigma_iX, converged = converged))
}
