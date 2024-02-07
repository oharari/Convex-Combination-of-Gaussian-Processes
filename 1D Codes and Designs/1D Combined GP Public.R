#**************************************************************************************
#***************************** 1D Combined GP code ************************************
#************************** Written by Ofir Harari ************************************
#*************** Uses Bayesian MCMC methods to fit a convex combination ***************
#************ of Gaussian Processes to the output of computer experiments *************
#**************************************************************************************




#**************************************************************************************
#************* Before running the script, make sure all of these packages *************
#************************ have been installed on your machine *************************
#**************************************************************************************
library(MASS)
library(mnormt)
library(LearnBayes)
library(coda)
library(CGP)
library(R.utils)
library(lhs)
library(MCMCpack)
library(ggplot2)
library(reshape2)
library(dplyr)
library(lemon)
library(latex2exp)
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
#******************* the set of univariate functions used in the paper ****************
#***************** the reader may add his own by typing "5" = ...... ******************
#**************************************************************************************
f <- function(x, code)
{
	switch(code,
	"1" = 0.5*sin(10*x)+0.5*cos(15*x),
	"2" = sin(10*x),
	"3" = sin(18*x-9)/(18*x-9),
	"4" = exp(3*x)*cos(5*(x-2))
	)
}
#**************************************************************************************



#**************************************************************************************
#**************************** Matern correlation function *****************************
#**** nu is the smothness parameter, theta the scale parameter and h the difference ***
#**************************************************************************************
Matern.corr.func <- function(nu, h, theta)
{
	ifelse(h==0,1,(2*sqrt(nu)*abs(h)/theta)^nu*besselK(2*sqrt(nu)*abs(h)/theta,nu)/(gamma(nu)*2^(nu-1)))
}
#**************************************************************************************



#**************************************************************************************
#********************* Vectorizing the Matern correlation function ********************
#**************************************************************************************
v.corr.func <- Vectorize(Matern.corr.func)
#**************************************************************************************



#**************************************************************************************
#**************** Evaluating the Gram matrix R for a given training set X *************
#****************** smoothness parameter nu and scale parameter theta *****************
#**************************************************************************************
corr.matrix <- function(nu, X, theta)
{
	n <- dim(X)[1]	
	A <- matrix(rep(X, n), nrow = n, byrow = T)
	U <- abs(A - t(A))
	matrix(v.corr.func(U, nu = nu, theta = theta), nrow = n)
}
#**************************************************************************************



#**************************************************************************************
#************* Evaluating the vector of correlations r between a new site x ***********
#********************************** and the training set X ****************************
#**************************************************************************************
corr.vec <- function(x, X, theta, nu)
{
	n <- dim(X)[1]
	A <- matrix(rep(x, n), nrow = n)
	U <- abs(A-X)
	v.corr.func(U, nu = nu, theta = theta)
}
#**************************************************************************************



#**************************************************************************************
#***** Maximum likelihood estimate for the intercept in the ordinary Kriging model ****
#**************************************************************************************
beta.MLE <- function(R.Inv, y)
{
	n <- length(c(y))
	t(as.vector(rep(1,n)))%*%R.Inv%*%c(y)/sum(R.Inv)
}
#**************************************************************************************




#**************************************************************************************
#**** Maximum likelihood estimate for the marginal variance in the ordinary Kriging ***
#********************* model (serves for the combined GP model too) *******************
#**************************************************************************************
sigma2.MLE <- function(R.Inv, y, beta)
{
	n <- length(y)
	t(y - c(beta)*rep(1,n))%*%R.Inv%*%(y - c(beta)*rep(1,n))/n
}
#**************************************************************************************




#**************************************************************************************
#********** The log-likelihood function, given R and the variance parameter ***********
#**************************************************************************************
log.like <- function(R, sigma2)
{
	n <- dim(R)[1]
	log(det(R, tol = 1e-16)) + n*log(sigma2)
}
#**************************************************************************************




#**************************************************************************************
#***************** The log-likelihood function, given all the parameters **************
#**************************************************************************************
log.likeli <- function(nu, theta, D.train, y.train)
{
	R <- corr.matrix(nu, D.train, theta)
	R.Inv <- solve(R, tol = 1e-16)
	beta <- beta.MLE(R.Inv, y.train)
	sigma2 <- sigma2.MLE(R.Inv, y.train, beta)
	return(log.like(R, sigma2))
}
#**************************************************************************************




#**************************************************************************************
#********* Finding the ML estimates by iteratively maximizing the ll function *********
#***** For given variance and intercept, using the resulting correlation parameter ****
#******************** theta for beta and sigma^2, and vice versa **********************
#**************************************************************************************
MLEs <- function(D, y, nu){
	theta.MLE <- NA
	while(is.na(theta.MLE[1])){	
		theta.MLE <- try(nlminb(runif(1) , log.likeli, nu = nu, D.train = D, 
		                        y.train = y)$par, 
		                 silent = TRUE)
		if(inherits(theta.MLE, "try-error")) theta.MLE <- NA
		else {
			R.Inv <- try(solve(corr.matrix(nu, D, theta.MLE)))
			if(inherits(R.Inv, "try-error")) theta.MLE <- NA
		}
	}
	beta.MLE <- beta.MLE(R.Inv, y)
	sigma2.MLE <- sigma2.MLE(R.Inv, y, beta.MLE)

	return(list(beta = beta.MLE, sigma2 = sigma2.MLE, theta = theta.MLE))
}
#**************************************************************************************




#**************************************************************************************
#* Calculating the posterior variance at the test set D.new for the ordinary Kriging **
#********* model based on the data collected at D.train and the ML estimates **********
#**************************************************************************************
post.var.single <- function(MLE, R.Inv, D.train, D.new, y.train, nu)
{
	n.train <- dim(R.Inv)[1]
	n.test <- dim(D.new)[1]

	theta <- MLE$theta
	sigma2 <- MLE$sigma2

	r <- apply(as.matrix(D.new), MARGIN = 1, FUN = corr.vec, X = D.train, theta = theta, nu = nu)
	var <- c(sigma2)*(rep(1, n.test) - diag(r%*%R.Inv%*%t(r)) + (rep(1, n.test) - c(t(rep(1, n.train))%*%R.Inv%*%t(r)))^2/sum(R.Inv))
	return(var)
}
#**************************************************************************************




#**************************************************************************************
#***** Calculating the posterior standard deviation at the test set D.new for the *****
#* ordinary Kriging model based on the data collected at D.train and the ML estimates *
#***** corrected for unknown sigma^2 (using unbiased estimation, as in Williams, ******
#***************************** Santner & Notz, 2003) **********************************
#**************************************************************************************
post.stdev.single <- function(MLE, R.Inv, D.train, D.new, y.train, nu)
{
	n.train <- dim(D.train)[1]
	sigma2 <- MLE$sigma2

	var <- post.var.single(MLE, R.Inv, D.train, D.new, y.train, nu)
	U <- apply(R.Inv, 1, sum)
	S <- matrix(rep(U, n.train), nrow = n.train, ncol = n.train, byrow = 1)/sum(U)
	Q.sq <- t(c(y.train))%*%(R.Inv -  R.Inv%*%S)%*%c(y.train)
	nu <- n.train - 1
	var.post <- c(Q.sq)*var/c(sigma2*nu)
	
	return(sqrt(var.post))
}
#**************************************************************************************




#**************************************************************************************
#******** Prediction Intervals for the ordinary Kriging model at new data D.new *******
#************************ Based on student's t distribution ***************************
#**************************************************************************************
CIs.single <- function(MLE, R.Inv, D.train, D.new, y.train, y.pred, alpha, nu)
{
	n.train <- dim(D.train)[1]
	stdev <- post.stdev.single(MLE, R.Inv, D.train, D.new, y.train, nu)
	df <- n.train - 1

	delta <- qt(1 - alpha/2, df)*stdev
	LL.Single <- y.pred - delta
	UL.Single <- y.pred + delta
	single.var <- stdev^2	

	return(data.frame(cbind(LL.Single, UL.Single, single.var)))
}
#**************************************************************************************




#**************************************************************************************
#******** Predictions and PIs for the ordinary Kriging model at new data D.new ********
#**************************************************************************************
prediction.single <- function(D.train, D.new, y.train, MLE, alpha, nu)
{	
	n.train <- dim(D.train)[1]
	n.new <- dim(D.new)[1]

	theta.MLE <- c(MLE$theta)
	beta.MLE <- as.numeric(MLE$beta)

	R.Inv <- solve(corr.matrix(nu, D.train, theta.MLE), solve = 1e-16)
	
	
	pred <- function(x)
	{
		beta.MLE + t(corr.vec(x, D.train, theta.MLE, nu=nu))%*%R.Inv%*%(y.train - beta.MLE*rep(1,n.train))
	}

	y.hat.single <- apply(D.new, 1, pred)
	CIs <- CIs.single(MLE, R.Inv, D.train, D.new, y.train, y.hat.single, alpha, nu)
	return(data.frame(cbind(y.hat.single, CIs)))
}
#**************************************************************************************




#**************************************************************************************
#**** Evaluating the Gram matrix R for a given training set X  for the Combined GP ****
#****************** model, given the parameters p, theta1 and theta2 ******************
#**************************************************************************************
Mixed.corr.matrix <- function(D.train, p, theta1, theta2, nu)
{
	R1 <- corr.matrix(nu, D.train, theta1)
	R2 <- corr.matrix(nu, D.train, theta2)
	R.mixed <- (p^2*R1 + (1-p)^2*R2)/(p^2 + (1-p)^2)

	return(R.mixed)
}
#**************************************************************************************



#**************************************************************************************
#*********** Evaluating the vector of correlations r for a new site x.new and *********
#********************************** taining set D.train *******************************
#**************************************************************************************
Mixed.corr.vec <- function(x.new, D.train, p, theta1, theta2, nu)
{
	corr1 <- corr.vec(x.new, D.train, theta1, nu=nu)
	corr2 <- corr.vec(x.new, D.train, theta2, nu=nu)

	return((p^2*corr1 + (1-p)^2*corr2)/(p^2 + (1-p)^2))
}
#**************************************************************************************




#**************************************************************************************
#*** The joint log-posterior density of p, theta1 and theta2 given the training data **
#**** The function returns the nominal value as well as the intercept and R matrix ****
#**************************************************************************************
logpost <- function(D.train, theta, y, sigma2, nu)
{
	psi1 <- theta[1]
	psi2 <- theta[2]
	phi <- theta[3]
	theta1 <- exp(psi1)
	theta2 <- exp(psi2)
	p <- 1/(1 + exp(-phi))

	R <- Mixed.corr.matrix(D.train, p, theta1, theta2, nu=nu)
		
	R.Inv <- NA

	R.Inv <- try(solve(R), silent = TRUE)
	if(inherits(R.Inv, "try-error")) R.Inv <- NA

	y <- c(y)
	beta <- beta.MLE(R.Inv, y)
	
	#*********** the multivariate normal log-likelihood *************
	log.like <- dmnorm(y, sum(beta), (p^2+(1-p)^2)*sigma2*R, log=1)

	#*** the log-jacobian of the transformation of the parameters ***
	#***** (done so they are supported by the entire real axis) *****
	log.jacob <- -phi -2*log(1 + exp(-phi)) + psi1 + psi2

      #******* the log-prior density (go ahead and change it!) ********
	log.prior <- - 4*psi1 - 2/theta1 - 6*psi2 - 16/theta2 #- 0.5*log(p*(1-p))

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
Metro <- function(start, N, samp.size, batch.size, alpha, D.train, sigma2, y, nu)
{
	n <- dim(D.train)[1]

	pb <- txtProgressBar(min = 0, max = N, style = 3)


	logpost.val <- function(theta)
	{
		return(logpost(D.train, theta, y, sigma2, nu)$val)
	}

	est <- laplace(logpost.val, start)
	pars <- list(mu = est$mode, v = est$var)

	samp <- matrix(0, ncol = length(start), nrow = N)

	beta <- list()
	R.Inv <- list()
	
	theta.old <- pars$mu


	k = 1
	l.old <- logpost(D.train, theta.old, y, sigma2, nu)
	pv = 0

	while(k <= N & pv < alpha)
	{		
		u <- runif(1)
		test <- matrix(NA, nrow = n, ncol = n)
		while(is.na(test[1,1]))
		{
			theta.candidate <- rmnorm(1, as.vector(theta.old), sqrt(2)*pars$v)		
			l.cand <- logpost(D.train, theta.candidate, y, sigma2, nu)
			test <- l.cand$R.Inv
		}

		R <- l.cand$val - l.old$val
		if(R > log(u))
		{
			samp[k,] <- theta.candidate
			theta.old <- theta.candidate
			beta <- c(beta, l.cand$beta)
			R.Inv <- c(R.Inv, list(test))
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

	return(list(sample=data.frame(samp[(k-samp.size):(k-1),]),
             beta=beta[(k-samp.size):(k-1)],R.Inv=R.Inv[(k-samp.size):(k-1)]))
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
factors.frame <- function(start, N, samp.size, batch.size, alpha.geweke, D.train, sigma2, y.train, net.samp.size, nu)
{
	n.train <- dim(D.train)[1]

	s <- Metro(start, N, samp.size, batch.size, alpha.geweke, D.train, sigma2, y.train, nu)
	temp <- data.frame(s$samp[(samp.size-net.samp.size+1):samp.size,])

	tempo <- temp[,1:3]
	names(tempo) <- c("log.theta1", "log.theta2", "logit.p")

	scoob <- mcmc(tempo)	

	dev.new(width = 18, height = 10)
	par(mfrow = c(2,3), mar = c(4,4,2,2))
      autocorr.plot(scoob[,1],auto.layout=F, lwd=2, main='Autocorrelation of '~log(theta[1]), cex.main=1.5)
	autocorr.plot(scoob[,2],auto.layout=F, lwd=2, main='Autocorrelation of '~log(theta[2]), cex.main=1.5)
	autocorr.plot(scoob[,3],auto.layout=F, lwd=2, main='Autocorrelation of '~logit(p), cex.main=1.5)

	traceplot(scoob[,1], main = 'Trace of '~log(theta[1]), cex.main=1.5)
	traceplot(scoob[,2], main = 'Trace of '~log(theta[2]), cex.main=1.5)
	traceplot(scoob[,3], main = 'Trace of '~logit(p), cex.main=1.5)
	
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
predict.post <- function(x.new, D.train, pars, sigma2, nu)
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

	r <- matrix(Mixed.corr.vec(x.new, D.train, p, theta1, theta2, nu=nu), ncol = 1)

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
prediction <- function(x.new, alpha, code, pars.frame, D.train, sigma2, nu)
{
	pred.samp <- data.frame(t(apply(pars.frame, 1, predict.post, x.new = x.new, 
                                D.train = D.train, sigma2 = sigma2, nu=nu)))
	n.samp <- dim(pred.samp)[1]
	names(pred.samp) <- c("mean", "var")

	y.hat.Combined <- mean(pred.samp$mean)
	y.true <- f(x.new, code = code)

	Var <- pred.samp$var
	stdev <- sqrt(Var)
	posterior.predictive <- rnorm(n.samp, pred.samp$mean, stdev)
	quants <- quantile(posterior.predictive, c(alpha/2, 1-alpha/2))
	Quant.Combined <- mean(y.true <= posterior.predictive)
	LL.Combined <- quants[1]
	UL.Combined <- quants[2]

	return(c(y.hat.Combined, Quant.Combined, LL.Combined, UL.Combined))
}
#**************************************************************************************




#**************************************************************************************
#******* Fitting both ordinary Kriging and combined GP and returning a table for ******
#********************************* comparison of the two ******************************
#**************************************************************************************
compare.GP <- function(D.new, alpha, code, params, D.train, sigma2, y.train, nu)
{
	n.train <- dim(D.train)[1]
	d <- dim(D.train)[2]
	n.new <- dim(D.new)[1]

	cgp <- CGP(X = D.train, yobs = c(y.train))

	pred.cgp <- predict(cgp, D.new, PI = TRUE)
	y.hat.CGP <- pred.cgp$Yp
	LL.CGP <- pred.cgp$Y_low
	UL.CGP <- pred.cgp$Y_up

	MLE <- MLEs(D.train, y.train, nu = nu)

	predict <- data.frame(t(apply_pb(D.new, 1, prediction, alpha = alpha, code = code, 
                              pars.frame = params, D.train = D.train, sigma2 = sigma2, nu=nu)))
	names(predict) <- c("y.hat.Combined", "Quant.Combined", "LL.Combined", "UL.Combined")
	y.true <- apply(D.new, 1, f, code = code)
	pred <- prediction.single(D.new = D.new, D.train = D.train, y.train = y.train, MLE = MLE, alpha = alpha, nu=nu)

	return(data.frame(cbind(X = D.new[,1], predict, pred, y.hat.CGP, LL.CGP, UL.CGP, y.true)))
}
#**************************************************************************************




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
	Coverage.Single <- mean((comp.object$y.true >= comp.object$LL.Single) & (comp.object$y.true <= comp.object$UL.Single))
	Coverage.CGP <- mean((comp.object$y.true >= comp.object$LL.CGP) & (comp.object$y.true <= comp.object$UL.CGP))

	Average.Quantile <- mean(comp.object$Quant.Combined)

	return(data.frame(cbind(RMSPE.Combined, RMSPE.Single, RMSPE.CGP, Coverage.Combined, 
             Coverage.Single, Coverage.CGP, Average.Quantile)))
}
#**************************************************************************************



prior_posterior_plot = function(x, theta1.prior, theta2.prior, params){
  thata1_lim = max(quantile(theta1.prior, .999),
                   quantile(params$theta1, 1))
  x1 = x[x <= thata1_lim]
  
  thata2_lim = max(quantile(theta2.prior, .999),
                   quantile(params$theta2, 1))
  x2 = x[x <= thata2_lim]
  
  df = params %>% 
    dplyr::select(p, theta1, theta2) %>% 
    dplyr::rename('$\\theta_1$' = 'theta1', 
                  '$\\theta_2$' = 'theta2', '$p$' = 'p') %>% 
    reshape2::melt() 
  
  l = lapply(as.list(unique(df$variable)), 
             function(x){
               d = density(df %>% 
                             filter(variable == x) %>% 
                             as.data.frame() %>% 
                             dplyr::select(value)%>% 
                             unlist()) 
               
               return(data.frame(Parameter = x, x = d$x, y = d$y))
             })
  
  df_posterior = do.call(rbind, l) %>% 
    mutate(Distribution = 'Posterior')
  
  p_lims = df_posterior %>% 
    filter(Parameter == '$p$') %>% 
    dplyr::select(x) %>% 
    range()
  
  df_prior = data.frame(Parameter = c(rep('$p$', 2), 
                                      rep('$\\theta_1$', length(x1)),
                                      rep('$\\theta_2$', length(x2))),
                        x = c(p_lims[1], p_lims[2], x1, x2),
                        y = c(1, 1, theta1.prior[1:length(x1)], 
                              theta2.prior[1:length(x2)]),
                        Distribution = 'Prior')
  
  df_all = rbind(df_prior, df_posterior) %>% 
    mutate(Distribution = factor(Distribution,
                                 levels = c('Prior', 'Posterior')))
  
  appender <- function(string){TeX(string)}
  
  p = ggplot(df_all, aes(x, y, col = Distribution, fill = Distribution)) + 
    geom_line() + 
    geom_area(stat = 'identity', alpha = .5) + 
    facet_rep_wrap(~ Parameter, scales = 'free', 
                   labeller = as_labeller(appender, 
                                          default = label_parsed)) + 
    theme(panel.background = element_blank(),
          legend.position = 'top',
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=16, hjust = .5,
                                    face="bold.italic"),
          strip.text.x = element_text(size = 12, face = "bold")) + 
    scale_y_continuous(expand = expansion(c(0, .025))) + 
    scale_x_continuous(expand = expansion(c(0, 0))) + 
    scale_color_brewer(palette = 'Set1') + 
    scale_fill_brewer(palette = 'Set1')
  
  return(p)
}




#**************************************************************************************
#******** This function receives training set 'D.train', y values 'y.train',  *********
#***** and returns a table of predictions and PIs for both ordinary Kriging and *******
#********************************* the Combined GP model ******************************
#**************************************************************************************
Combined.GP.fit <- function(D.train, y.train, D.new, start, N.max, samp.size, 
                            alpha.geweke, batch.size, nu, code)
{
	test <- data.frame(matrix(NA, nrow = dim(D.new)[1], ncol = 4))
	while(is.na(test[1,1])){
		MLE <- MLEs(D.train, y.train, nu=nu)
		sigma2 <- as.numeric(MLE$sigma2)
			
		test <- try(prediction.single(D.train, D.new, c(y.train), MLE, alpha, nu), 
		            silent = TRUE)
		if(inherits(test, "try-error")) 
		  test <- data.frame(matrix(NA, nrow = dim(D.new)[1], ncol = 4))
	}	

	print("MCMC Progress:", quote = 0)
	params <- factors.frame(start, N.max, samp.size, batch.size, alpha.geweke, 
	                        D.train, sigma2, y.train, net.samp.size, nu)

	x = seq(0.001, 40, by = 0.005)
	theta1.prior = dinvgamma(x, 3, 2)
	theta2.prior = dinvgamma(x, 5, 16)
	p = prior_posterior_plot(x, theta1.prior, theta2.prior, params)
  print(p)

	print("Prediction Progress:", quote = 0)
	comp.object <- compare.GP(D.new, alpha, code, params, D.train, sigma2, y.train, nu)

	return(list(comp.object = comp.object, D.train = D.train, y.train = y.train))
			
}
#**************************************************************************************




#**************************************************************************************
#******** Plotting the resulting Combined GP and ordinary Kriging predictions *********
#**************************************************************************************
plot.GP <- function(Data, F_Num)
{
	comp.object <- Data$comp.object
	D.train <- Data$D.train
	y.train <- Data$y.train
	code <- F_Num
	l <- dim(comp.object)[1]	

	x.sim <- seq(0, 1, length = 100)
	true.func <- f(x.sim, code)

	ylims <- c(
 			min(c(min(true.func),min(unlist(comp.object$LL.Combined)))), 
            	max(c(max(true.func), max(unlist(comp.object$UL.Combined))))
                 )

	dev.new(width = 16, height = 8)
	par(mar = c(4,5,1,1))
	plot(x.sim, true.func, type = 'l', lwd = 2, col = 2, xlab = '', ylab = '', axes = 0,
           ylim = ylims, xlim = c(0,1))
	polygon(c(comp.object$X, sort(comp.object$X, decreasing = T)), 
              c(comp.object$LL.Combined, comp.object$UL.Combined[l:1]), col = 'light gray', border = NA)
	lines(x.sim, true.func, lwd = 2, col = 2)
	lines(comp.object$X, comp.object$y.hat.single, col = 4, lwd = 5, cex = 1.5, lty = 3)
	lines(comp.object$X, comp.object$y.hat.Combined, lty = 2, lwd = 3)
	axis(1, at = seq(-1, 1, by = .2), cex.axis=1.3)
	axis(2, at = round(seq(min(comp.object$LL.Combined)-1, max(comp.object$UL.Combined)+1, by  = 0.3), 1), cex.axis=1.3)
	mtext(expression(x), 1, 2, cex = 1.8)
	mtext(expression(y(x)), 2, 3, cex = 1.8)
	points(D.train, y.train, pch = 19, cex = 2.5, col = 'dark green')

	if(F_Num==1|F_Num==2)
	{
		legend('bottomleft', c("True", "Combined", "Single"), lty = c(1,2,3), col = c(2,1,4), lwd = c(3,3,5), cex = 1.7)
	}
	if(F_Num==3|F_Num==4)
	{
		legend('topleft', c("True", "Combined", "Single"), lty = c(1,2,3), col = c(2,1,4), lwd = c(3,3,5), cex = 1.7)
	}
}
#**************************************************************************************





#**************************************************************************************
#******************************** End of Functions!!! *********************************
#***************************** Simulation starts here!!! ******************************
#**************************************************************************************


n.train <- 8 #training set size
nu <- 5 #smoothness parameter for the Matern correlation function
d <- 1 #dimensionality of the problem
alpha <- 0.05 #complement of the coverage probability of the prediction intervals
start <- c(0,1.5,0) # a starting point for the Markov chain for the MCMC sampling
N.max <- 10000 #maximum number of iterations before terminating the MCMC process
samp.size <- 5000 #chain length to perform the Geweke test for convergence to the stationary distribution
net.samp.size <- 2500 #sample size from the posterior to retain
alpha.geweke <- 0.5 #confidence level for the Geweke test
batch.size <- 20 #number of samples between Geweke tests for convergence to the stationary


n.new <- 50 #test set size
D.new <- as.matrix(seq(0, 1, length = n.new), ncol = 1) #test set

func <- 1 #function number (originally 1-4)
D.train <- randomLHS(n.train, 1) #training set
y.train <- f(D.train, func) #training 'y' values
data <- Combined.GP.fit(D.train, y.train, D.new, start, N.max, samp.size, 
                        alpha.geweke, batch.size, nu, func)
Comparison.Summary(data$comp.object)
plot.GP(data, func)


