#**************************************************************************************
#********* 2D Combined GP code with Batch-Sequential Maximum Entropy Designs **********
#************************** Written by Ofir Harari ************************************
#*************** Uses Bayesian MCMC methods to fit a convex combination ***************
#******** of Isotropic Gaussian Processes to the output of computer experiments *******
#**************************************************************************************




#**************************************************************************************
#************* Before running the script, make sure all of these packages *************
#************************ have been installed on your machine *************************
#**************************************************************************************
library(MASS)
library(mnormt)
library(LearnBayes)
library(coda)
library(R.utils)
library(lhs)
library(MCMCpack)
library(mlegp)
library(CGP)
library(rgl)
library(pscl)
#**************************************************************************************




#*************************************************************************************
#********* a technical function used to print a progress bar on the screen ***********
#*************************************************************************************
apply_pb <- function(X, MARGIN, FUN, ...)
{
	env <- environment()
	pb_Total <- sum(dim(X)[MARGIN])
	counter <- 0
	pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

	wrapper <- function(...)
	{
		curVal <- get("counter", envir = env)
    		assign("counter", curVal +1 ,envir= env)
		setTxtProgressBar(get("pb", envir= env), curVal +1)
		FUN(...)
	}
	res <- apply(X, MARGIN, wrapper, ...)
	close(pb)
	res
}
#**************************************************************************************




#*************************************************************************************
#*********** a bug-free version of the CGP function from the CGP package *************
#*************************************************************************************
CGP = function (X, yobs, nugget_l = 0.001, num_starts = 5, theta_l = NULL, 
                alpha_l = NULL, kappa_u = NULL) 
{
  yobs <- as.numeric(yobs)
  DD <- as.matrix(X)
  var_names <- colnames(DD)
  n <- nrow(DD)
  p <- ncol(DD)
  one <- rep(1, n)
  onem <- rep(1, n - 1)
  Stand_DD <- apply(DD, 2, function(x) (x - min(x))/max(x - 
                                                          min(x)))
  scales <- apply(DD, 2, function(x) max(x) - min(x))
  if (length(theta_l) > 1) 
    print("Error: Lower bound for theta needs to be specified as a scalar!")
  if (length(alpha_l) > 1) 
    print("Error: Lower bound for alpha needs to be specified as a scalar!")
  if (length(kappa_u) > 1) 
    print("Error: Upper bound for kappa needs to be specified as a scalar!")
  if (is.null(theta_l)) 
    theta_l <- 1e-04
  theta_lower <- rep(theta_l, p)
  if (is.null(alpha_l)) 
    alpha_l <- log(10^2) * mean(1/dist(Stand_DD)^2)
  theta_upper <- rep(alpha_l, p)
  kappa_l <- alpha_l
  if (is.null(kappa_u)) 
    kappa_u <- log(10^6) * mean(1/dist(Stand_DD)^2)
  if (sum(theta_l > alpha_l) > 0) 
    print("Error: Lower bound of theta exceeds the upper bound!")
  lower = c(nugget_l, theta_lower, kappa_l, 0)
  upper = c(1, theta_upper, kappa_u, 1)
  PSI <- function(theta) {
    A <- DD %*% diag(sqrt(theta), ncol = p)
    A <- as.matrix(dist(A, diag = T, upper = T))
    R <- exp(-A^2)
    return(R)
  }
  Stand_PSI <- function(Stand_theta) {
    A <- Stand_DD %*% diag(sqrt(Stand_theta), ncol = p)
    A <- as.matrix(dist(A, diag = T, upper = T))
    R <- exp(-A^2)
    return(R)
  }
  var.MLE.DK <- function(ww) {
    lambda <- ww[1]
    Stand_theta <- ww[2:(p + 1)]
    kappa <- ww[p + 2]
    Stand_alpha <- kappa + Stand_theta
    bw <- ww[p + 3]
    G <- Stand_PSI(Stand_theta)
    L <- Stand_PSI(Stand_alpha)
    Gbw <- Stand_PSI(Stand_theta * bw)
    Sig <- diag(n)
    for (rep in 1:4) {
      Q <- G + lambda * Sig^(1/2) %*% L %*% Sig^(1/2)
      invQ <- solve(Q)
      beta <- c((one %*% invQ %*% yobs)/(one %*% invQ %*% 
                                           one))
      temp <- invQ %*% (yobs - beta * one)
      gip <- beta * one + G %*% temp
      e <- yobs - gip
      Sig <- diag(c(Gbw %*% e^2/(Gbw %*% one)))
      Sig2 <- mean(diag(Sig))
      Sig <- Sig/Sig2
    }
    Q <- G + lambda * Sig^(1/2) %*% L %*% Sig^(1/2)
    invQ <- solve(Q)
    beta <- c((one %*% invQ %*% yobs)/(one %*% invQ %*% one))
    tau2 <- t(yobs - beta * one) %*% invQ %*% (yobs - beta * 
                                                 one)/n
    val <- c(log(det(Q)) + n * log(tau2))
    if (!is.finite(val)) 
      val <- 1e+06
    return(val)
  }
  
  n_par <- p + 3
  n_candidate <- 500 + num_starts
  LHD <- function(N, k) {
    x <- matrix(rep(1:N, k), ncol = k, nrow = N)
    x <- apply(x, 2, sample)
    x <- (x - 0.5)/N
    return(x)
  }
  
  starts <- LHD(n_candidate, n_par)
  range_S <- matrix(rep(upper - lower, n_candidate), ncol = n_par, 
                    byrow = TRUE)
  low_S <- matrix(rep(lower, n_candidate), ncol = n_par, byrow = TRUE)
  starts <- starts * range_S + low_S
  cand_obj <- apply(starts, 1, var.MLE.DK)
  index <- (rank(cand_obj, ties.method = "min") <= num_starts)
  Starts <- starts[index, ]
  Fop <- function(x) optim(x, var.MLE.DK, lower = lower, upper = upper, 
                           method = "L-BFGS-B")$val
  op_obj <- apply(Starts, 1, Fop)
  beststart <- Starts[which.min(op_obj), ]
  op <- optim(beststart, var.MLE.DK, lower = lower, upper = upper, 
              method = "L-BFGS-B")
  objval <- op$val
  lambda <- op$par[1]
  Stand_theta <- op$par[2:(p + 1)]
  kappa <- op$par[p + 2]
  Stand_alpha <- kappa + Stand_theta
  bw <- op$par[p + 3]
  theta <- Stand_theta/scales^2
  alpha <- Stand_alpha/scales^2
  Yjfp <- rep(0, n)
  for (jf in 1:n) {
    G <- PSI(theta)[-jf, -jf]
    L <- PSI(alpha)[-jf, -jf]
    Gbw <- PSI(theta * bw)[-jf, -jf]
    Sig <- diag(n - 1)
    for (rep in 1:4) {
      Q <- G + lambda * Sig^(1/2) %*% L %*% Sig^(1/2)
      invQ <- solve(Q)
      beta <- c((onem %*% invQ %*% yobs[-jf])/(onem %*% invQ %*% 
                                                 onem))
      temp <- invQ %*% (yobs[-jf] - beta * onem)
      gip <- beta * onem + G %*% temp
      e <- yobs[-jf] - gip
      Sig <- diag(c(Gbw %*% e^2/(Gbw %*% onem)))
      Sig2 <- mean(diag(Sig))
      Sig <- Sig/Sig2
    }
    Q <- G + lambda * Sig^(1/2) %*% L %*% Sig^(1/2)
    invQ <- solve(Q)
    beta <- c((onem %*% invQ %*% yobs[-jf])/(onem %*% invQ %*% 
                                               onem))
    tau2 <- t(yobs[-jf] - beta * onem) %*% invQ %*% (yobs[-jf] - 
                                                       beta * onem)/(n - 1)
    temp <- invQ %*% (yobs[-jf] - beta * onem)
    g <- PSI(theta)[jf, -jf]
    l <- PSI(alpha)[jf, -jf]
    gbw <- PSI(theta * bw)[jf, -jf]
    vjf <- (t(gbw) %*% (e^2)/(t(gbw) %*% onem))/Sig2
    vjf <- as.vector(vjf)
    q <- g + lambda * sqrt(vjf) * Sig^(1/2) %*% l
    Yjfp[jf] <- beta + t(q) %*% temp
  }
  rmscv <- sqrt(sum((yobs - Yjfp)^2)/n)
  G <- PSI(theta)
  L <- PSI(alpha)
  Gbw <- PSI(theta * bw)
  Sig <- diag(n)
  for (rep in 1:4) {
    Q <- G + lambda * Sig^(1/2) %*% L %*% Sig^(1/2)
    invQ <- solve(Q)
    beta <- c((one %*% invQ %*% yobs)/(one %*% invQ %*% one))
    temp <- invQ %*% (yobs - beta * one)
    gip <- beta * one + G %*% temp
    e <- yobs - gip
    res2 <- e^2
    Sig <- diag(c(Gbw %*% res2/(Gbw %*% one)))
    sf <- mean(diag(Sig))
    Sig <- Sig/sf
  }
  Q <- G + lambda * Sig^(1/2) %*% L %*% Sig^(1/2)
  invQ <- solve(Q)
  beta <- c((one %*% invQ %*% yobs)/(one %*% invQ %*% one))
  tau2 <- t(yobs - beta * one) %*% invQ %*% (yobs - beta * 
                                               one)/n
  temp <- invQ %*% (yobs - beta * one)
  theta <- matrix(theta, nrow = 1)
  alpha <- matrix(alpha, nrow = 1)
  if (!is.null(var_names)) {
    colnames(theta) = var_names
    colnames(alpha) = var_names
  }
  est <- list(X = DD, yobs = yobs, var_names = var_names, lambda = lambda, 
              theta = theta, alpha = alpha, bandwidth = bw, Sig_matrix = Sig, 
              sf = sf, res2 = res2, temp_matrix = temp, invQ = invQ, 
              mu = beta, tau2 = tau2, beststart = beststart, objval = objval, 
              rmscv = rmscv, Yp_jackknife = Yjfp)
  est$call <- match.call()
  class(est) <- "CGP"
  return(est)
}
#*************************************************************************************




#*************************************************************************************
#******* a bug-free version of the predict.CGP function from the CGP package *********
#*************************************************************************************
predict.CGP <- function(object, newdata=NULL, PI=FALSE,...){
  
  UU<-newdata
  DD<-object$X
  yobs<-object$yobs
  n<-nrow(DD)
  p<-ncol(DD)
  one<-rep(1,n)
  
  lambda<-object$lambda
  theta<-as.vector(object$theta)
  alpha<-as.vector(object$alpha)
  bw<-object$bandwidth
  Sig<-object$Sig_matrix
  sf<-object$sf
  res2<-object$res2
  temp<-object$temp_matrix
  invQ<-object$invQ
  tau2<-c(object$tau2)
  beta<-object$mu
  
  if(is.null(UU)){
    Yp<-NULL
    gp<-NULL
    lp<-NULL
    v<-NULL
    Y_low<-NULL
    Y_up<-NULL
  }
  
  if(!is.null(UU)){
    UU<-as.matrix(UU)
    N<-nrow(UU)
    if(ncol(UU)!=p) print("Predictive location input UU is of wrong dimension!")
    ppp<-rep(0,N)
    g<-rep(0,n)
    gbw<-rep(0,n)
    l<-rep(0,n)
    Yp<-rep(0,N)
    gp<-rep(0,N)
    lp<-rep(0,N)
    v<-rep(0,N)
    for (k in 1:N){
      for (r in 1:n){
        g[r]<-exp(-(DD[r,]-UU[k,])^2%*%(theta))
        gbw[r]<-exp(-(DD[r,]-UU[k,])^2%*%(theta*bw))
        l[r]<-exp(-(DD[r,]-UU[k,])^2%*%(alpha))
      }
      v[k]<-(t(gbw)%*%(res2)/(t(gbw)%*%one))/sf
      q<-g+lambda*sqrt(v[k])*Sig^(1/2)%*%l
      Yp[k]<-beta+t(q)%*%temp
      gp[k]<-beta+t(g)%*%temp
      if(PI){
        lp[k]<-lambda*sqrt(v[k])*t(l)%*%Sig^(1/2)%*%temp
        ppp[k]<-1+lambda*v[k]-t(q)%*%invQ%*%q+(1-t(q)%*%invQ%*%one)^2/(one%*%invQ%*%one)
      }
    }
    if(PI){
      ppp[ppp<0]<-0
      ka<-1.96
      Y_up<-Yp+ka*sqrt(tau2*ppp)
      Y_low<-Yp-ka*sqrt(tau2*ppp)
    }
    if(!PI){
      Y_low<-NULL
      Y_up<-NULL
    }
  }
  
  val<-list(Yp=Yp,gp=gp,lp=lp,v=v,Y_low=Y_low,Y_up=Y_up)
  return(val)
  
}
#*************************************************************************************




#**************************************************************************************
#******************* the set of bivariate functions used in the paper *****************
#***************** the reader may add his own by typing "6" = ...... ******************
#**************************************************************************************
f <- function(vars, code)
{
	x <- vars[1]
	y <- vars[2]
	switch(code,
		 "1" = ((x-0.2)^2 - (y-0.7)^2)*exp(-5*((x-0.8)^2 + (y-0.1)^2))*cos(10*(x-0.5)*y),
		 "2" = ((x-0.5)^2 + 4*(y-0.8)^2)*(cos(pi*(x-0.1)) + cos(pi*(y-0.5)))
	)
}
#**************************************************************************************




#**************************************************************************************
#**************** Evaluating the Gram matrix R for a given training set X *************
#******************* and scale parameter theta (Isotropic process) ********************
#**************************************************************************************
corr.matrix.ISO <- function(X, theta)
{
	n <- dim(X)[1]
	d <- dim(X)[2]
	Theta <- diag(rep(theta, d))
	U <- matrix(apply(X^2%*%Theta, 1, sum), nrow = n, ncol = n ,byrow=F)
	V <- -2*X%*%Theta%*%t(X)
	Dist <- U + t(U) + V
	return(exp(-Dist))
}
#**************************************************************************************




#**************************************************************************************
#************* Evaluating the vector of correlations r between a new site x ***********
#******************** and the training set X (Isotropic process) **********************
#**************************************************************************************
corr.vec.ISO <- function(x, X, theta)
{
	n <- dim(X)[1]
	d <- dim(X)[2]
	theta <- rep(theta, d)
	Theta <- diag(theta)
	r <- as.vector(exp(-(rep(t(theta)%*%x^2, n) - 2*X%*%Theta%*%x + apply(X^2%*%Theta, 1, sum))))
	return(r)	
}
#**************************************************************************************




#**************************************************************************************
#***** Maximum likelihood estimate for the intercept in the ordinary Kriging model ****
#**************************************************************************************
beta.MLE <- function(R.Inv, y)
{
	n <- length(y)
	c(t(rep(1,n))%*%R.Inv%*%y/sum(R.Inv))
}
#**************************************************************************************




#**************************************************************************************
#**** Evaluating the Gram matrix R for a given training set X  for the Combined GP ****
#****************** model, given the parameters p, theta1 and theta2 ******************
#**************************************************************************************
Mixed.corr.matrix <- function(D.train, p, theta1, theta2)
{
	R1 <- corr.matrix.ISO(D.train, theta1)
	R2 <- corr.matrix.ISO(D.train, theta2)
	R.mixed <- (p^2*R1 + (1-p)^2*R2)/(p^2 + (1-p)^2)

	return(R.mixed)
}
#**************************************************************************************




#**************************************************************************************
#*********** Evaluating the vector of correlations r for a new site x.new and *********
#********************************** taining set D.train *******************************
#**************************************************************************************
Mixed.corr.vec <- function(x.new, D.train, p, theta1, theta2)
{
	corr1 <- corr.vec.ISO(x.new, D.train, theta1)
	corr2 <- corr.vec.ISO(x.new, D.train, theta2)

	return((p^2*corr1 + (1-p)^2*corr2)/(p^2 + (1-p)^2))
}
#**************************************************************************************




#**************************************************************************************
#*** The joint log-posterior density of p, theta1 and theta2 given the training data **
#**** The function returns the nominal value as well as the intercept and R matrix ****
#**************************************************************************************
logpost <- function(D.train, theta, y, sigma2)
{
	psi1 <- theta[1]
	psi2 <- theta[2]
	phi <- theta[3]
	theta1 <- exp(psi1)
	theta2 <- exp(psi2)
	p <- 1/(1 + exp(-phi))

	R <- Mixed.corr.matrix(D.train, p, theta1, theta2)
		
	R.Inv <- NA

	R.Inv <- try(solve(R), silent = TRUE)
	if(inherits(R.Inv, "try-error")) R.Inv <- NA

	beta <- beta.MLE(R.Inv, y)

	log.like <- dmnorm(y, sum(beta), (p^2+(1-p)^2)*sigma2*R, log=1)
	log.jacob <- -phi -2*log(1 + exp(-phi)) + psi1 + psi2
	log.prior <- - 4*psi1 - 2/theta1 - 6*psi2 - 16/theta2 # - 0.5*log(p*(1-p))
	val <- log.like + log.jacob + log.prior

	return(list(val=val, beta=beta, R.Inv = R.Inv))
}
#**************************************************************************************





#**************************************************************************************
#************** This is the heart of the script - the metropolis sampling *************
#************** from the posterior distribution of p, theta1 and theta2 ***************
#*** The proposal density is multivariate normal based on the laplace approximation ***
#**** of the posterior, and after 'batch.size' samples are gathered, a Geweke test ****
#*** is performed to check for stationarity of the Markov chain, until 'samp.size' ****
#*** samples are gathered. The procedure halts after a maximum of N iterations, and ***
#**** the last 'samp.size' samples are retained. Along with the random sample, the ****
#** already inverted R matrices and beta intercepts are returned (for each triplet) *** 
#**************************************************************************************
Metro <- function(start, N, samp.size, batch.size, alpha, D.train, sigma2, y, ...)
{
	pb <- txtProgressBar(min = 0, max = N, style = 3)


	logpost.val <- function(theta)
	{
		return(logpost(D.train, theta, y, sigma2)$val)
	}

	est <- laplace(logpost.val, start)
	pars <- list(mu = est$mode, v = est$var)

	samp <- matrix(0, ncol = length(start), nrow = N)

	beta <- list()
	R.Inv <- list()
	log.posterior <- list()
	
	theta.old <- pars$mu


	k = 1
	l.old <- logpost(D.train, theta.old, y, sigma2, ...)
	pv = 0

	while(k <= N & pv < alpha)
	{		
		u <- runif(1)
		theta.candidate <- rmnorm(1, as.vector(theta.old), sqrt(2)*pars$v)
		l.cand <- logpost(D.train, theta.candidate, y, sigma2, ...)
		R <- l.cand$val - l.old$val

		if(R > log(u))
		{
			samp[k,] <- theta.candidate
			theta.old <- theta.candidate
			beta <- c(beta, l.cand$beta)
			R.Inv <- c(R.Inv, list(l.cand$R.Inv))
			log.posterior <- c(log.posterior, l.cand$val)
			
			l.old <- l.cand
			k <- k+1
			setTxtProgressBar(pb, k)
		}

		if((k-1) >= samp.size & (k-1)%%batch.size == 0)
		{
			options(warn=2)
			pv = try(min(2*(1-pnorm(abs(geweke.diag(mcmc(samp[(k-samp.size):(k-1)]))$z)))),silent=TRUE)
			if(inherits(pv, "try-error")) pv <- 0
		}				
	}

	close(pb)

	return(list(sample = data.frame(samp[(k-samp.size):(k-1),]),
             beta = beta[(k-samp.size):(k-1)], R.Inv = R.Inv[(k-samp.size):(k-1)],
             logpost = log.posterior[(k-samp.size):(k-1)]))
}
#**************************************************************************************




#**************************************************************************************
#****** pre-calculating some of the terms that will later be used for prediction ******
#****************************** for the combined GP model *****************************
#**************************************************************************************
factors <- function(MCMC.data, n.train, y.train)
{
	R.Inv <- matrix(as.numeric(MCMC.data[1:n.train^2]), nrow = n.train)
	beta <- as.numeric(MCMC.data[n.train^2+1])
	mean.factor <- as.vector(R.Inv%*%(y.train - rep(beta, n.train)))
	var.factor1 <- apply(R.Inv, 2, sum)
	var.factor2 <- sum(var.factor1)

	return(c(mean.factor, var.factor1, var.factor2))
}
#**************************************************************************************




#**************************************************************************************
#****** Given a starting point (for the Markov chain to roll on), training data *******
#** and the Matern smoothness parameter, this function returns a random sample from ***
#********** the posterior distribution, along with the pre-calculated terms ***********
#******************************** for future predictions ******************************
#**************************************************************************************
factors.frame <- function(start, N, samp.size, batch.size, alpha.geweke, D.train, sigma2, y.train, net.samp.size)
{
	n.train <- dim(D.train)[1]

	s <- Metro(start, N, samp.size, batch.size, alpha.geweke, D.train, sigma2, y.train)
	temp <- data.frame(s$samp[(samp.size-net.samp.size+1):samp.size,])
	p <- 1/(1 + exp(-temp[,3]))
	theta1 <- exp(temp[,1])
	theta2 <- exp(temp[,2])
	samp <- data.frame(cbind(p, theta1, theta2))
	names(samp) <- c("p", "theta1", "theta2")
	beta <- as.numeric(unlist(s$beta))[(samp.size-net.samp.size+1):samp.size]
	names(beta) <- "beta"
	R.Inv <- data.frame(t(matrix(unlist(s$R.Inv), ncol = samp.size))[(samp.size-net.samp.size+1):samp.size,])

	prediction.factors <- data.frame(t(apply(data.frame(cbind(R.Inv, beta)), 1, factors, n.train = n.train, y.train = y.train)))
	return(cbind(samp, beta, prediction.factors, R.Inv))
}
#**************************************************************************************





#**************************************************************************************
#**** Given a new site x.new, training data and a random triplet from the posterior ***
#******** this function returns a pair consisting of the mean and variance of *********
#***************** the posterior predictive distribution of y at x.new ****************
#**************************************************************************************
predict.post <- function(x.new, D.train, pars, sigma2)
{
	n <- dim(D.train)[1]
	p <- as.numeric(pars[1])
	theta1 <- as.numeric(pars[2])
	theta2 <- as.numeric(pars[3])
	beta <- as.numeric(pars[4])
	mean.factor <- as.numeric(pars[5:(4+n)])
	var.factor1 <- as.numeric(pars[(5+n):(4+2*n)])
	var.factor2 <- as.numeric(pars[(5+2*n)])
	R.Inv <- matrix(as.numeric(pars[(6+2*n):(5+2*n+n^2)]), nrow=n)	

	r <- matrix(Mixed.corr.vec(x.new, D.train, p, theta1, theta2), ncol = 1)

	var <- as.numeric(sigma2*(1-t(r)%*%R.Inv%*%r+(1-t(var.factor1)%*%r)^2/var.factor2))
	mean <- as.numeric(beta + t(as.vector(mean.factor))%*%r)	
		
	return(cbind(mean, var))
}
#**************************************************************************************





#**************************************************************************************
#**** This function uses the random sample drawn from the posterior distribution of ***
#***** p, theta1 and theta2, to draw a random sample from the posterior predictive ****
#***** distribution of y at x.new, and returns the average and suitable percentiles ***
#***** for prediction and prediction intervals, as well as the average quantile  ******
#*** of the true y value within the posterior predictive distribution of y at x.new ***
#**************************************************************************************
prediction <- function(x.new, alpha, code, pars.frame, D.train, sigma2)
{
	pred.samp <- data.frame(t(apply(pars.frame, 1, predict.post, x.new = x.new, D.train = D.train, sigma2 = sigma2)))
	names(pred.samp) <- c("mean", "var")
	pred.samp <- pred.samp[pred.samp$var>=0, ]
	n.samp <- dim(pred.samp)[1]

	mu0 <- mean(pred.samp$mean)
	Var <- pred.samp$var
	stdev <- sqrt(Var)
	posterior.predictive <- rnorm(n.samp, pred.samp$mean, stdev)
	quants <- quantile(posterior.predictive, c(alpha/2, 1-alpha/2))
	Quant1 <- mean(mu0 <= posterior.predictive)
	LL1 <- quants[1]
	UL1 <- quants[2]	


	return(c(mu0, Quant1, LL1, UL1))
}
#**************************************************************************************




#**************************************************************************************
#********* Fitting ordinary Kriging, combineg GP and CGP (Ba & Joseph 2012) ***********
#***************** and returning a table for comparison of the three ******************
#**************************************************************************************
compare.GP <- function(D.new, alpha, code, params, D.train, sigma2, y.train)
{
	n.train <- dim(D.train)[1]
	d <- dim(D.train)[2]
	n.new <- dim(D.new)[1]

	cgp <- CGP(D.train, y.train)
	ord.krig <- mlegp(D.train, y.train)

	pred <- predict.CGP(cgp, D.new, PI = TRUE)
	y.hat.CGP <- pred$Yp
	LL.CGP <- pred$Y_low
	UL.CGP <- pred$Y_up

	pred.single <- predict.gp(ord.krig, D.new, se.fit = TRUE)
	sigma <- pred.single$se.fit
	y.hat.single <- pred.single$fit
	LL.single <- pred.single$fit - sigma*qt(1-alpha/2, df=n.train-1)
	UL.single <- pred.single$fit + sigma*qt(1-alpha/2, df=n.train-1)

	predict <- data.frame(t(apply_pb(D.new, 1, prediction, alpha = alpha, code = code, 
                              pars.frame = params, D.train = D.train, sigma2 = sigma2)))
	names(predict) <- c("y.hat.Combined", "Quant.Combined", "LL.Combined", "UL.Combined")
	y.true <- apply(D.new, 1, f, code = code)

	return(data.frame(cbind(X1 = D.new[,1], X2 = D.new[,2], predict, y.hat.single, LL.single, 
                        UL.single, y.hat.CGP, LL.CGP, UL.CGP, y.true)))
}


#**************************************************************************************
#********** Summary Statistics for the comparison of the ordinary Kriging *************
#**** and Combined GP models, given a table produced by the 'compare.GP' function *****
#**************************************************************************************
Comparison.Summary <- function(comp.object)
{
	e.Combined <- comp.object$y.true - comp.object$y.hat.Combined
	e.Single <- comp.object$y.true - comp.object$y.hat.single
	e.CGP <- comp.object$y.true - comp.object$y.hat.CGP

	RMSPE.Combined <- sqrt(mean(e.Combined^2))
	RMSPE.Single <- sqrt(mean(e.Single^2))
	RMSPE.CGP <- sqrt(mean(e.CGP^2))
	
	Coverage.Combined <- mean((comp.object$y.true >= comp.object$LL.Combined) & (comp.object$y.true <= comp.object$UL.Combined))
	Coverage.Single <- mean((comp.object$y.true >= comp.object$LL.single) & (comp.object$y.true <= comp.object$UL.single))
	Coverage.CGP <- mean((comp.object$y.true >= comp.object$LL.CGP) & (comp.object$y.true <= comp.object$UL.CGP))

	Average.Quantile <- mean(comp.object$Quant.Combined)

	return(data.frame(cbind(RMSPE.Combined, RMSPE.Single, RMSPE.CGP, Coverage.Combined, 
             Coverage.Single, Coverage.CGP, Average.Quantile)))
}
#**************************************************************************************




#**************************************************************************************
#********** Summary Statistics for the comparison of the ordinary Kriging *************
#**** and Combined GP models, given a table produced by the 'compare.GP' function *****
#**************************************************************************************
Comparison.Summary <- function(comp.object)
{
	comp.object <- comp.object[!is.na(comp.object$y.true),]
	e.Combined <- comp.object$y.true - comp.object$y.hat.Combined
	e.Single <- comp.object$y.true - comp.object$y.hat.single
	e.CGP <- comp.object$y.true - comp.object$y.hat.CGP

	RMSPE.Combined <- sqrt(mean(e.Combined^2))
	RMSPE.Single <- sqrt(mean(e.Single^2))
	RMSPE.CGP <- sqrt(mean(e.CGP^2))
	
	Coverage.Combined <- mean((comp.object$y.true >= comp.object$LL.Combined) & (comp.object$y.true <= comp.object$UL.Combined))
	Coverage.Single <- mean((comp.object$y.true >= comp.object$LL.single) & (comp.object$y.true <= comp.object$UL.single))
	Coverage.CGP <- mean((comp.object$y.true >= comp.object$LL.CGP) & (comp.object$y.true <= comp.object$UL.CGP))

	Average.Quantile <- mean(comp.object$Quant.Combined)

	return(data.frame(cbind(RMSPE.Combined, RMSPE.Single, RMSPE.CGP, Coverage.Combined, 
             Coverage.Single, Coverage.CGP, Average.Quantile)))
}
#**************************************************************************************




#**************************************************************************************
#************* Fitting the Combined GP model and making predictions *******************
#**************************************************************************************
Combined.GP.fit <- function(D.train, y.train, D.new, start, N.max, samp.size, 
                            alpha.geweke, batch.size, code)
{
	ord <- mlegp(D.train, y.train)		
	sigma2 <- ord$sig2

	print("MCMC Progress:", quote = 0)
	params <- factors.frame(start, N.max, samp.size, batch.size, alpha.geweke, D.train, 
                              sigma2, y.train, net.samp.size)

	theta1 <- density(params[,2])
	theta2 <- density(params[,3])
	x <- seq(0.01, 20, by = 0.01)
	theta1.prior <- dinvgamma(x, 3, 2)
	theta2.prior <- dinvgamma(x, 5, 16)
	p <- density(params[,1])

	ylims1 <- c(0, max(c(theta1$y, theta1.prior)))
	ylims2 <- c(0, max(c(theta2$y, theta2.prior)))

	dev.new(width = 15, height = 5)
	par(mfrow = c(1,3))

	plot(p, lwd = 2, main = '')
	title(main = expression(paste("Distribution of ", p)), cex.main = 1.5)
	lines(c(0, 1), c(1, 1), col = 2, lty = 2, lwd = 2)
	legend('topright', col=c(2,1), lty = c(2,1), c("Prior", "Posterior"), lwd = 2, 
              cex = 1.2)
			
	plot(theta1, lwd = 2, main = '', xlim = c(0,6), ylim = ylims1)
	title(main = expression(paste("Distribution of ", theta[1])), cex.main = 1.5)
	lines(x, theta1.prior, lwd = 2, col = 2, lty = 2)
	legend('topright', col=c(2,1), lty = c(2,1), c("Prior", "Posterior"), lwd = 2, 
              cex = 1.2)

	plot(theta2, lwd = 2, main = '', ylim = ylims2, xlim = c(0,25))
	title(main = expression(paste("Distribution of ", theta[2])), cex.main = 1.5)
	lines(x, theta2.prior, lwd = 2, col = 2, lty = 2)
	legend('topright', col=c(2,1), lty = c(2,1), c("Prior", "Posterior"), lwd = 2, 
             cex = 1.2)


	print("Prediction Progress:", quote = 0)
	comp.object <- compare.GP(D.new, alpha, code, params, D.train, sigma2, y.train)

	return(comp.object)
			
}
#**************************************************************************************




#**************************************************************************************
#****************************** Plotting the resulting fit ****************************
#**************************************************************************************
plot2dGP <- function(Data, f_num)
{
	Data.to.show <- Data
	n <- dim(Data.to.show)[1]
	Data.to.show <- Data.to.show[sample(1:n, size = 120, replace = 0),]

	x1 <- Data.to.show$X1
	x2 <- Data.to.show$X2
	mu <- Data.to.show$y.hat.Combined

	x <- y <- seq(-1, 1, length = 50)
	f2 <- function(x, y, code) 
	{
		switch(code,
			"1" = ((x-0.2)^2 - (y-0.7)^2)*exp(-5*((x-0.8)^2 + (y-0.1)^2))*cos(10*(x-0.5)*y),
			"2" = ((x-0.5)^2 + 4*(y-0.8)^2)*(cos(pi*(x-0.1)) + cos(pi*(y-0.5)))
			)
	}

	z <- outer(x, y, f2, code = paste(f_num))

	Interval.data <- data.frame(rbind(cbind(x1, x2, Data.to.show$LL.Combined), cbind(x1, x2, Data.to.show$UL.Combined)))
	Interval.data <- Interval.data[order(Interval.data$x1, Interval.data$x2),]

	open3d(windowRect=c(100,100,1000,750))
	bg3d("white")
	material3d(col="black")	
	persp3d(x, y, z, aspect=c(1, 1, 1), col = "lightgreen", axes = 0, xlab = "", ylab = "", zlab = "", smooth = T, type = "s")
	points3d(x1, x2, mu, size = 6)
	segments3d(Interval.data, col = 'red')	
}
#**************************************************************************************




#**************************************************************************************
#***************** Cross correlation matrix between old and new inputs ****************
#**************************************************************************************
cross.corr.matrix <- function(D.old, D.new, theta)
{
	n1 <- dim(D.new)[1]
	n2 <- dim(D.old)[1]
	d <- dim(D.new)[2]
	Theta <- diag(rep(theta,d))
	
	U <- matrix(apply(D.new^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=F)
	V <- -2*D.new%*%Theta%*%t(D.old)
	W <- matrix(apply(D.old^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=T)
	Dist <- U + V + W

	return(exp(-Dist))
}
#**************************************************************************************



#**************************************************************************************
#********************* The Entropy criterion for the first batch **********************
#**************************************************************************************
Entropy <- function(D, p, theta1, theta2)
{
	R <- Mixed.corr.matrix(D, p, theta1, theta2)

	return(-det(R))
}
#**************************************************************************************



#**************************************************************************************
#******************* The Entropy criterion for two batches of inputs ******************
#**************************************************************************************
Augmented.Mixed.Entropy <- function(D.old, D.new, p, theta1, theta2, R.old.Inv)
{
	n.new <- dim(D.new)[1]
	R.cross <- (p^2*cross.corr.matrix(D.old, D.new, theta1) + (1-p)^2*cross.corr.matrix(D.old, D.new, theta2))/(p^2 + (1-p)^2)
	R.new <- Mixed.corr.matrix(D.new, p, theta1, theta2)
	E <- try(-det(R.new - R.cross%*%R.old.Inv%*%t(R.cross), tol = 1e-16), silent = T)

	return(E)
}
#**************************************************************************************




#**************************************************************************************
#**************** Finding the batch-sequential Maximum Entropy Design *****************
#**************************************************************************************
Entropy.optim <- function(n, d, p, theta1, theta2, n.starts)
{
	vals <- rep(0, n.starts)
	designs <- list()

	Entropy.univar <- function(D)
	{
		D <- matrix(D, ncol = d)
		Entropy(D, p, theta1, theta2)
	}
	
	for(k in 1:n.starts)
	{
		start <- c(-1+2*optimumLHS(n,d))
		optimum <- optim(start, Entropy.univar, method = "L-BFGS-B", lower = rep(-1,n*d), 
                             upper = rep(1,n*d))
		D.optimum <- matrix(optimum$par, ncol = d)
		Entropy.optimum <- optimum$value
		designs[[k]] <- D.optimum
		vals[k] <- Entropy.optimum
	}

	min.val <- vals[which.min(vals)]
	best.design <- designs[[which.min(vals)]]
	
	return(list(Design = best.design, log.entropy = -min.val))
}
#**************************************************************************************



#**************************************************************************************
#***************** Finding the batch-sequentialMaximum Entropy Design *****************
#**************************************************************************************
Batch.Entropy.optim <- function(D.old, n.new, d, p, theta1, theta2, n.starts)
{
	vals <- rep(0, n.starts)
	designs <- list()
	R.old <- (p^2*corr.matrix.ISO(D.old, theta1) + (1-p)^2*corr.matrix.ISO(D.old, theta2))/(p^2 + (1-p)^2)
	R.old.Inv <- solve(R.old, tol = 1e-16)

	Entropy.univar <- function(D.new)
	{
		D.new <- matrix(D.new, ncol = d)
		Augmented.Mixed.Entropy(D.old, D.new, p, theta1, theta2, R.old.Inv)
	}
	
	for(k in 1:n.starts)
	{
		start <- c(-1+2*optimumLHS(n.new,d))
		optimum <- optim(start, Entropy.univar, method = "L-BFGS-B", lower = rep(-1,n.new*d), 
                             upper = rep(1,n.new*d))
		D.optimum <- matrix(optimum$par, ncol = d)
		Entropy.optimum <- optimum$value
		designs[[k]] <- D.optimum
		vals[k] <- Entropy.optimum
	}

	min.val <- vals[which.min(vals)]
	best.design <- designs[[which.min(vals)]]
	
	return(list(Design = best.design, log.entropy = -min.val))
}
#**************************************************************************************







#**************************************************************************************
#******************************** End of Functions!!! *********************************
#***************************** Simulation starts here!!! ******************************
#**************************************************************************************


#**************************************************************************************
#************* path to the first batch (size 14) Maximum Entropy design ***************
wd = dirname(rstudioapi::getActiveDocumentContext()$path)
#**************************************************************************************
setwd(wd)
#**************************************************************************************


#**************************************************************************************
func <- 1 #function number
alpha <- 0.05 #complement of the coverage probability of the prediction intervals
start <- c(0,.5,1) #a starting point for the Markov chain for the MCMC sampling
N.max <- 5000 #maximum number of iterations before terminating the MCMC process
samp.size <- 1000 #chain length to perform the Geweke test for convergence to the stationary distribution
net.samp.size <- 1000 #sample size from the posterior to retain
alpha.geweke <- 0.5 #confidence level for the Geweke test
batch.size <- 20 #number of samples between Geweke tests for convergence to the stationary
n.old <- 14 #first batch size
p.prior <- 0.5
theta1.prior <- 1 #qigamma(.5, 3, 2)
theta2.prior <- 4#qigamma(.5, 5, 16)
n.starts <- 20 #number of starting positions for the optimization stage
#Calculating the first batch design, based on the prior medians
#D.old <- Entropy.optim(n.old, d=2, p.prior, theta1.prior, theta2.prior, n.starts)$Design
D.old <- as.matrix(read.table("Initial ME Design.txt"))
plot(D.old, cex = 3, pch = 19, xlim = c(-1,1), ylim = c(-1,1), main = "First Batch Design",
     cex.main = 2, xlab = expression(X[1]), ylab = expression(X[2]), cex.lab = 1.6)

#**************************************************************************************


y.old <- apply(D.old, 1, f, code = paste(func))
ord <- mlegp(D.old, y.old)		
sigma2 <- ord$sig2
n.new <- 7 #second batch size

#**************************************************************************************
#**************** Sampling from the posterior to derive MAP estimators ****************
#**************************************************************************************
params <- Metro(start, N.max, 2000, batch.size, alpha.geweke, D.old, sigma2, y.old)
p.old <- 1/(1 + exp(-params$sample[,3]))
theta1.old <- exp(params$sample[,1])
theta2.old <- exp(params$sample[,2])
pars.old <- cbind(p.old, theta1.old, theta2.old)
log.post <- unlist(params$logpost)
pars.MAP <- pars.old[which.max(log.post),]

#p <- pars.MAP[1]
#theta1 <- pars.MAP[2]
#theta2 <- pars.MAP[3]

p <- mean(p.old)
theta1 <- mean(theta1.old)
theta2 <- mean(theta2.old)
#**************************************************************************************


#**************************************************************************************
#************** Calculating 2nd batch design based on the MAP estimators **************
#**************************************************************************************
D.optim <- Batch.Entropy.optim(D.old, n.new, 2, p, theta1, theta2, 25)
R.old <- (p^2*corr.matrix.ISO(D.old, theta1) + (1-p)^2*corr.matrix.ISO(D.old, theta2))/(p^2 + (1-p)^2)

val <- round(D.optim$log.entropy + log(det(R.old)), 3)
dev.new(width = 8, height = 6)
par(mar = c(4,0,4,0))
plot(D.old, pch = 19, cex = 2, col = 1, xlab = "", ylab = "", 
     main = bquote(paste("Maximum Entropy Design, criterion value = ",.(val))),
     xlim = c(-1.9,1.9), ylim = c(-1.1,1.1), axes = 0, cex.main = 1.6)
lines(c(-1.1,-1.1,1.1,1.1,-1.1), c(-1.1,1.1,1.1,-1.1,-1.1), lwd = 2)
points(D.optim$Design, pch = 21, cex = 2.5, lwd = 3, bg = 'light gray')
mtext(bquote(paste(hat(p)==.(round(p,3)), " , ", hat(theta)[1]==.(round(theta1,3)), " , ", hat(theta)[2]==.(round(theta2,3)))),
      1, 1, cex = 1.5)
legend(1.1, 0.3, c("1st Batch", "2nd Batch"), pch = c(19, 21), pt.lwd= 3, pt.bg = c('black', 'light gray'),
       pt.cex = 2, cex= 1.3)
#**************************************************************************************


#**************************************************************************************
#****************** Sampling from the posterior after the 2nd batch *******************
#**************************************************************************************
D.train <- as.matrix(rbind(D.old, D.optim$Design))
y.train <- apply(D.train, 1, f, code = paste(func))
ord <- mlegp(D.train, y.train)		
sigma2 <- ord$sig2
params <- Metro(start, N.max, samp.size, batch.size, alpha.geweke, D.train, sigma2, y.train)
p.new <- 1/(1 + exp(-params$sample[,3]))
theta1.new <- exp(params$sample[,1])
theta2.new <- exp(params$sample[,2])
#**************************************************************************************


#**************************************************************************************
#********************** Plotting the posteriors after each stage **********************
#**************************************************************************************
p.old <- density(p.old)
p.new <- density(p.new)
theta1.prior <- density(1/rgamma(100000, 3, 2))
theta1.old <- density(theta1.old)
theta1.new <- density(theta1.new)
theta2.prior <- density(1/rgamma(100000, 5, 16))
theta2.old <- density(theta2.old)
theta2.new <- density(theta2.new)
ylim1 <- max(c(p.old$y, p.new$y))
ylim2 <- max(c(theta1.prior$y, theta1.old$y, theta1.new$y))
ylim3 <- max(c(theta2.prior$y, theta2.old$y, theta2.new$y))

dev.new(width = 18, height = 6)
par(mfrow = c(1,3))
plot(c(0,1), c(1,1), main = "", lwd = 3, ylim = c(0,ylim1), type = 'l', xlab = '', ylab = '',
     cex.axis = 1.3)
lines(p.old, col = 2, lty = 2, lwd =3)
lines(p.new, col = 4, lty = 3, lwd = 5)
title(main = expression(paste(scriptstyle(p),"|", scriptstyle(y), sep = "")), cex.main = 3)
legend('topright', c("Prior", "1st Posterior", "2nd Posterior"), col = c(1,2,4), lty = c(1,2,3),
       lwd = c(3,3,4), cex = 1.5)

plot(theta1.prior, lwd = 2, xlim = c(0,8), main = "", ylim = c(0,ylim2), xlab = '', ylab = '',
     cex.axis = 1.3)
lines(theta1.old, col = 2, lwd = 3, lty = 2)
lines(theta1.new, col = 4, lty = 3, lwd = 5)
title(main = expression(paste(scriptstyle(theta[1]),"|", scriptstyle(y), sep = "")), cex.main = 3)
legend('topright', c("Prior", "1st Posterior", "2nd Posterior"), col = c(1,2,4), lty = c(1,2,3),
       lwd = c(3,3,4), cex = 1.5)

plot(theta2.prior, lwd = 2, xlim = c(0,20), main = "", sub = "", ylim = c(0,ylim3), xlab = '', 
     ylab = '', cex.axis = 1.3)
lines(theta2.old, col = 2, lwd = 3, lty  = 2)
lines(theta2.new, col = 4, lty = 3, lwd = 5)
title(main = expression(paste(scriptstyle(theta[2]),"|", scriptstyle(y), sep = "")), cex.main = 3)
legend('topright', c("Prior", "1st Posterior", "2nd Posterior"), col = c(1,2,4), lty = c(1,2,3),
       lwd = c(3,3,4), cex = 1.5)
#**************************************************************************************


#**************************************************************************************
#***************** Fitting the Combined GP model on the entire design *****************
#**************************************************************************************
n.new <- 25 #test set size (per axis)
u <- seq(-1, 1, length = n.new)
D.new <- expand.grid(u,u) #test set (nXn lattice)
data <- Combined.GP.fit(D.train, y.train, D.new, start, N.max, samp.size, alpha.geweke,
                        batch.size, func)
Comparison.Summary(data)
plot2dGP(data, func)
#**************************************************************************************

