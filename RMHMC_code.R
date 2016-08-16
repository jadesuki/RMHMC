rm(list=ls())
#Attaching dependencies
library(SoDA)
library(MASS)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(inline)


#setting working directory 
setwd("~/Documents/Git/RMHMC")
#load simulated data 
load('sim_AUG.Rdata')
#loading functions 
source('RMHMC_function.R')


#C++ code for generating samples from high dimensional multivariate normal distribution
code <- '
using namespace Rcpp;
int n = as<int>(n_);
arma::vec mu = as<arma::vec>(mu_);
arma::mat sigma = as<arma::mat>(sigma_);
int ncols = sigma.n_cols;
arma::mat Y = arma::randn(n, ncols);
return wrap(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
'
#Creating function 
rmvnorm.rcpp <- 
  cxxfunction(signature(n_="integer", mu_="numeric", sigma_="matrix"), code,
              plugin="RcppArmadillo", verbose=TRUE)


#Setting data 
spiketrain = spiketrain[1:1000]
stimulus = stimulation[1:1000]


#Initialization of Markov Chain 
niter = 1000

beta.fix = 1
sigma.fix = 0.1

theta.mat = matrix(rep(0, 3 * niter), nrow = 3, ncol = niter) #ordering: mu, phi (gamma), alpha
state.mat = matrix(rep(0,length(stimulus)+1),nrow = length(stimulus)+1,ncol=niter)
theta.mat[,1] = c(-5, -1, 1.5)
state.mat[,1] = E.x(theta = theta.mat[,1], i.data = stimulus)


#Simulate Markov chain with RMHMC with Gibbs scheme
for(iter in 2:niter){
	print(iter)

	#Setting current theta and state 
	current.theta = theta.mat[,(iter-1)]
	current.state = state.mat[,(iter-1)]

	#Part 1: Update latent states
	#Getting tensor of states 
	tensor.state = tensorX(theta = current.theta, statevec = current.state)
	#Inverse of the tensor 
	inv.tensor.state = chol2inv(chol(tensor.state))
  	
  	#Generate auxiliary particles 
	current.aux = rmvnorm.rcpp(1,rep(0,length(current.state)),tensor.state)

	#running leapfrog method to generate candidates 
	leapfrog.X = gen.X(x.data = current.state, n.data=spiketrain, i.data = stimulus, theta= current.theta,aux = current.aux, Xten = tensor.state, Xten.inv=inv.tensor.state)
	(candidate.X = leapfrog.X$Xcandidate)
	candidate.p = leapfrog.X$aux 

	#Finding new tensor regarding candidates of latent states 
	new.tensor.state = tensorX(theta = current.theta, statevec = candidate.X)
	#Finding inverse of the new tensor 
	new.inv.tensort.state = chol2inv(chol(new.tensor.state))
	
	#Evaluate Hamiltonian function based on current states and  candidate states 
	logtemp1 = Ham.x(x.data = candidate.X,  n.data = spiketrain, i.data = stimulus, theta = current.theta, aux = candidate.p, Xten= new.tensor.state, Xten.inv = new.inv.tensort.state )
	logtemp2 = Ham.x(x.data = current.theta,n.data = spiketrain, i.data = stimulus, theta = current.theta, aux = current.aux, Xten=tensor.state,   Xten.inv = inv.tensor.state)

	logtemp3 = logtemp2 - logtemp1
	#Doing metropolis update 
	if(!is.na(logtemp3) && runif(1) < exp(logtemp3)){
		current.state = candidate.X; state.mat[,iter] = candidate.X
	}else{
	  state.mat[,iter] = current.state
	}

	#Part 2: Update theta 
	#Calculate mean vector 
	mean.states   = E.x(theta = current.theta, i.data = stimulus)
	#Calculate variance vector 
	var.states    = V.x(theta = current.theta, n = length(stimulus))

	#getting tenosr and inverse tensor of theta 
	tensor.theta = tensorPar(theta = current.theta, meanvec = mean.states, varvec = var.states, i.data = stimulus)
	inv.tensor.theta = solve(tensor.theta)

	#Generate auxiliary variables 
	current.aux = rmvnorm(1,rep(0,3),tensor.theta)

	#Simulate a trajectory from Hamiltonian system by leapfrog method 
	leapfrog.Par = gen.Par(x.data = current.state, n.data=spiketrain, i.data = stimulus, theta = current.theta, aux = current.aux ,meanvec = mean.states, varvec=var.states,Parten = tensor.theta, Parten.inv = inv.tensor.theta)
	#Getting candidates of theta and auxiliary 
	(candidate.Par = leapfrog.Par$Parcandidate)
	candidate.aux   = leapfrog.Par$aux 

	#Getting new mean vector and variance vector 
	new.mean.states = E.x(theta = candidate.Par, i.data = stimulus)
	new.var.states  = V.x(theta = candidate.Par, n = length(stimulus))

	#Getting new tenwor and its inverse 
	new.tensor.theta = tensorPar(theta = candidate.Par, meanvec  = new.mean.states, varvec = new.var.states,i.data = stimulus)
  new.inv.tensor.theta = solve(new.tensor.theta)
  
  	#Evaluate Hamiltonian function based on current theta and candidate theta 
	logtemp4 = Ham.par(x.data = current.state, n.data = spiketrain, i.data = stimulus, theta = candidate.Par,aux = candidate.aux,Parten = new.tensor.theta,Parten.inv = new.inv.tensor.theta)
	logtemp5 = Ham.par(x.data = current.state, n.data = spiketrain, i.data = stimulus, theta = current.theta,aux = current.aux,  Parten = tensor.theta, Parten.inv = inv.tensor.theta)

   
	logtemp6 = logtemp5 - logtemp4

	#Do Metropolis update 
	if(!is.na(logtemp6)&&runif(1)<exp(logtemp6)){
		current.theta = candidate.Par; theta.mat[,iter] = candidate.Par
	}else{
	  theta.mat[,iter] = current.theta
	}
}

