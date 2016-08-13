rm(list=ls())
library(SoDA)
library(MASS)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(inline)

setwd("~/Dropbox/RMHMC")
load('sim_AUG.Rdata')
source('RMHMC_function.R')

code <- '
using namespace Rcpp;
int n = as<int>(n_);
arma::vec mu = as<arma::vec>(mu_);
arma::mat sigma = as<arma::mat>(sigma_);
int ncols = sigma.n_cols;
arma::mat Y = arma::randn(n, ncols);
return wrap(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
'

rmvnorm.rcpp <- 
  cxxfunction(signature(n_="integer", mu_="numeric", sigma_="matrix"), code,
              plugin="RcppArmadillo", verbose=TRUE)


spiketrain = spiketrain[1:100]
stimulus = stimulation[1:100]
niter = 10000

beta.fix = 1
sigma.fix = 0.1


theta.mat = matrix(rep(0, 3 * niter), nrow = 3, ncol = niter) #ordering: mu, phi (gamma), alpha
state.mat = matrix(rep(0,length(stimulus)+1),nrow = length(stimulus)+1,ncol=niter)
theta.mat[,1] = c(-1, -1.8, 2)
state.mat[,1] = E.x(theta = theta.mat[,1], i.data = stimulus)
iter = 2
for(iter in 2:niter){
	print(iter)

	current.theta = theta.mat[,(iter-1)]
	current.state = state.mat[,(iter-1)]

	tensor.state = tensorX(theta = current.theta, statevec = latents[1:length(state.mat[,1])])
	inv.tensor.state = chol2inv(chol(tensor.state))
  
	current.aux = rmvnorm.rcpp(1,rep(0,length(current.state)),tensor.state)

	leapfrog.X = gen.X(x.data = current.state, n.data=spiketrain, i.data = stimulus, theta= current.theta,aux = current.aux, Xten = tensor.state, Xten.inv=inv.tensor.state)
	(candidate.X = leapfrog.X$Xcandidate)
	candidate.p = leapfrog.X$aux 

	new.tensor.state = tensorX(theta = current.theta, statevec = candidate.X)
	new.inv.tensort.state = chol2inv(chol(new.tensor.state))
	
	logtemp1 = Ham.x(x.data = candidate.X,  n.data = spiketrain, i.data = stimulus, theta = current.theta, aux = candidate.p, Xten= new.tensor.state, Xten.inv = new.inv.tensort.state )
	logtemp2 = Ham.x(x.data = current.theta,n.data = spiketrain, i.data = stimulus, theta = current.theta, aux = current.aux, Xten=tensor.state,   Xten.inv = inv.tensor.state)

	logtemp3 = logtemp2 - logtemp1

	if(!is.na(logtemp3) && runif(1) < exp(logtemp3)){
		current.state = candidate.X; state.mat[,iter] = candidate.X
	}else{
	  state.mat[,iter] = current.state
	}


	mean.states   = E.x(theta = current.theta, i.data = stimulus)
	var.states    = V.x(theta = current.theta, n = length(stimulus))
	tensor.theta = tensorPar(theta = current.theta, meanvec = mean.states, varvec = var.states, i.data = stimulus)
	inv.tensor.theta = solve(tensor.theta)

	current.aux = rmvnorm(1,rep(0,3),tensor.theta)

	leapfrog.Par = gen.Par(x.data = current.state, n.data=spiketrain, i.data = stimulus, theta = current.theta, aux = current.aux ,meanvec = mean.states, varvec=var.states,Parten = tensor.theta, Parten.inv = inv.tensor.theta)
	(candidate.Par = leapfrog.Par$Parcandidate)
	candidate.aux   = leapfrog.Par$aux 

	new.mean.states = E.x(theta = candidate.Par, i.data = stimulus)
	new.var.states  = V.x(theta = candidate.Par, n = length(stimulus))

	new.tensor.theta = tensorPar(theta = candidate.Par, meanvec  = new.mean.states, varvec = new.var.states,i.data = stimulus)
    new.inv.tensor.theta = solve(new.tensor.theta)
  
	logtemp4 = Ham.par(x.data = current.state, n.data = spiketrain, i.data = stimulus, theta = candidate.Par,aux = candidate.aux,Parten = new.tensor.theta,Parten.inv = new.inv.tensor.theta)
	logtemp5 = Ham.par(x.data = current.state, n.data = spiketrain, i.data = stimulus, theta = current.theta,aux = current.aux,  Parten = tensor.theta, Parten.inv = inv.tensor.theta)

    #logtemp4 = HMC.ham.Par(x.data = current.state, n.data = spiketrain, i.data = stimulus, theta = candidate.Par,aux = candidate.aux)
	#logtemp5 = HMC.ham.Par(x.data = current.state, n.data = spiketrain, i.data = stimulus, theta = current.theta,aux = current.aux)

	logtemp6 = logtemp5 - logtemp4

	if(!is.na(logtemp6)&&runif(1)<exp(logtemp6)){
		current.theta = candidate.Par; theta.mat[,iter] = candidate.Par
	}else{
	  theta.mat[,iter] = current.theta
	}
}

