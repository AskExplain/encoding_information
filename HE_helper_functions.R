## Kun Yue, yuek@uw.edu
## 02/24/2020

## main function for estimating linear mixed models via REHE
## Code is adapted from R package GENESIS (v2.14.3)

## For faster computation, set VConly to TRUE if only computing point estimations for variance components (OLS estimation for beta will be computed)
## If subsequent GWAS analysis will be performed, need to set VConly=F to obtain the required intermediate quantities (WLS estimation for beta will be computed)

## computeCI=T if would like to compute confidence intervals for variance components. Fixed effects will use weighted lease square.

fitREHE = function(y, X, K, group.idx=NULL, computeCI = F, VConly = F){ 
  
  # y is the response vector, X is the fixed covariate design matrix (the first column being 1), K is the list of correlation matrices
  
  
  if (is.null(group.idx)){group.idx <- list(resid.var = 1:length(y))}else{stop("Heterogeneous residual variance feature is not developed.")}
  
  XTXXT = solve(t(X)%*%X, t(X))
  beta = XTXXT%*%y
  # P = X%*%XTXXT
  n = length(y)
  In = diag(1, n)
  nK = length(K)
  
  residual = y - X%*%beta
  
  y_matrix = as.vector(tcrossprod(residual))
  X_matrix = cbind(as.vector(In), sapply(1:nK, function(w)as.vector(K[[w]])))
  
  XTX = crossprod(X_matrix)
  XTY = crossprod(X_matrix, y_matrix)
  
  
  res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,nK+1), bvec = rep(0, nK+1), meq=0, factorized=FALSE)
  
  estimates = res_qp$solution
  estimates_he = res_qp$unconstrained.solution
  
  
  # numerical correction for REHE estimations (due to machine error)
  for(s in 1:length(estimates)){
    estimates[s] = ifelse(estimates[s]<0, 0, estimates[s])}
  
  estimates = c(estimates[2:(nK+1)], residual=estimates[1])
  names(estimates) = c(paste('K', 1:nK), 'residual')
  
  if(!VConly || computeCI){ 
    # functions to prepare compatible quantities for GWAS analysis under GENESIS package
    sq = .computeSigmaQuantities(varComp = estimates, covMatList = K, group.idx = group.idx)
    lq = .calcLikelihoodQuantities(Y = y, X = X, Sigma.inv = sq$Sigma.inv, cholSigma.diag = sq$cholSigma.diag)
    output = list(varComp = estimates, AI = NA, converged = T, zeroFLAG = 0, niter = 0,
                  Sigma.inv = sq$Sigma.inv, beta = lq$beta, residM = lq$residM, fits = lq$fits, eta = 0, 
                  logLikR = lq$logLikR, logLik = lq$logLik, RSS = lq$RSS, In=In, cholSigma = sq$cholSigma)
  }else{ 
    sq=list(NULL)
    lq=list(NULL)
    output = list(varComp = estimates, beta = beta)
  }
  
  
  return(output)
  
  
}






.calcLikelihoodQuantities <- function(Y, X, Sigma.inv, cholSigma.diag){
  
  if (is(Sigma.inv, "Matrix")) X <- Matrix(X)
  n <- length(Y)
  k <- ncol(X)
  
  ### Calulate the weighted least squares estimate
  Sigma.inv_X <- crossprod(Sigma.inv, X)
  Xt_Sigma.inv_X <- crossprod(X, Sigma.inv_X)
  # fix issue with not recognizing the matrix as symmetric
  Xt_Sigma.inv_X <- (Xt_Sigma.inv_X + t(Xt_Sigma.inv_X))/2
  chol.Xt_Sigma.inv_X <- chol(Xt_Sigma.inv_X)
  Xt_Sigma.inv_X.inv <- chol2inv(chol.Xt_Sigma.inv_X)
  beta <- crossprod(Xt_Sigma.inv_X.inv, crossprod(Sigma.inv_X, Y))
  
  # calc Xb
  fits <- X %*% beta
  # calc marginal residuals = (Y - Xb)
  residM <- as.vector(Y - fits)
  
  ### calculate PY
  PY <- crossprod(Sigma.inv, residM)
  
  # compute RSS
  YPY <- crossprod(Y, PY)
  RSS <- as.numeric(YPY/(n-k))
  # Sigma.inv_R <- crossprod(Sigma.inv, residM)
  # Rt_Sigma.inv_R <- crossprod(residM, Sigma.inv_R)
  # Rt_Sigma.inv_R <- crossprod(residM, PY)
  # RSS <- as.numeric(Rt_Sigma.inv_R/(n - k)) 
  
  # log likelihood
  logLik <- as.numeric(-0.5 * n * log(2 * pi * RSS) - sum(log(cholSigma.diag)) - 0.5 * YPY/RSS)
  # REML log likelihood; accounting for estimation of mean effects 
  logLikR <- as.numeric(logLik + 0.5 * k * log(2 * pi * RSS) - sum(log(diag(chol.Xt_Sigma.inv_X))))
  
  return(list(PY = PY, RSS = RSS, logLik = logLik, logLikR = logLikR, 
              Sigma.inv_X = Sigma.inv_X, Xt_Sigma.inv_X.inv = Xt_Sigma.inv_X.inv, 
              beta = as.numeric(beta), fits = fits, residM = residM))
  
}

.computeSigmaQuantities <- function(varComp, covMatList, group.idx = NULL, vmu = NULL, gmuinv = NULL){
  m <- length(covMatList)
  n <- nrow(covMatList[[1]])
  
  # contribution to Sigma from random effects
  Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
  
  if (is.null(vmu)){
    # gaussian family
    # contribution to Sigma from residual variance
    if (is.null(group.idx)){
      diagV <- rep(varComp[m+1],n)
    } else{
      g <- length(group.idx)
      diagV <- rep(NA, n)
      for(i in 1:g){
        diagV[group.idx[[i]]] <- varComp[m+i]
      }
      # mylevels <- rep(NA, n)
      # for(i in 1:g){
      #     mylevels[group.idx[[i]]] <- i # setting up vector of indicators; who is in which group
      # }
      # diagV <- (varComp[(m+1):(m+g)])[mylevels]
    }
    
  } else {
    # non-gaussian family
    diagV <- as.vector(vmu)/as.vector(gmuinv)^2
  }
  
  # construct Sigma
  Sigma <- Vre
  diag(Sigma) <- diag(Sigma) + diagV   
  # cholesky decomposition
  cholSigma <- chol(Sigma)
  # inverse
  Sigma.inv <- chol2inv(cholSigma)
  
  return(list(Sigma.inv = Sigma.inv, Vre = Vre, diagV = diagV, cholSigma.diag = diag(cholSigma)))
  
}
