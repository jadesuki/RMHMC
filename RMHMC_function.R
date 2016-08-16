#Joint likelihood of latent states x_{0:k} and spiketrain N_{1:k} given theta 
like = function(log = TRUE, x.data, n.data , i.data ,theta, beta = beta.fix,sigma = sigma.fix,delta = 1) {

  mu = theta[1]
  phi = tanh(theta[2])
  alpha = theta[3]
  
  k = length(i.data)

  #seperate evaluation of joint likelihood 
  temp1 = -(k + 1) * log(2 * pi) / 2 - (k + 1) * log(sigma^2)

  temp2 = -sum((x.data[-1] - phi * x.data[-length(x.data)] - alpha * i.data)^2) / (2 * sigma^2)

  temp3 = log(1 - phi^2) / 2 - x.data[1]*(1 - phi^2) / (2 * sigma^2)

  temp4 = sum(n.data * (mu + beta * x.data[-1] + log(delta)) - exp(mu + beta * x.data[-1]) * delta)
  
  result = temp1 + temp2 + temp3 + temp4
  
  #returns log of joint likelihood 
  if (log) {return(result)}
  #returns joint likelihood
  else {return(exp(result))}
}


#Prior
#flat priors of latent statex 
#prior.x = function(log = TRUE, x.data){
#  priors = sum(dnorm(x.data,mean = 0 , sd = 10, log = TRUE))
#  if (log) {
#  	return(priors)
#  }else{
#  	return(exp(priors))
#  }
#}

#flat priors of theta 
prior.par = function(log = TRUE, theta, sd.prior = 10){
	priors = sum(dnorm(theta, mean = 0 , sd = sd.prior, log = TRUE))
	if(log){
		return(priors)
	}else{
		return(exp(priors))
	}
}



#Posterior
#Posterior of latent states is the joint likelihood of states and spiketrain 
post.x = function(log = TRUE, x.data,n.data,i.data, theta) {
  result = like(log = TRUE, x.data = x.data, n.data = n.data, i.data = i.data, theta = theta, delta = 1)

  if (log) {return(result)}
  else {return(exp(result))}
}

#Posterior of theta 
#with flat priors N(0,sd=10)
post.theta = function(log=TRUE, x.data,n.data, i.data,theta, beta=beta.fix, sigma=sigma.fix) {
	#log(posterior) = log(likelihood) + log(prior)
	results = like(log=TRUE, x.data = x.data,n.data = n.data,i.data = i.data, theta = theta,delta = 1) + prior.par(log = TRUE, theta = theta)

	if(log){
		return(results)
	}else{
		return(exp(results))
	}
}



#Expectation and Variance of Latent states xk's
#Getting mean vector of xk , E(xk)
E.x = function(theta, i.data){
  mu = theta[1]
  phi = tanh(theta[2])
  alpha = theta[3]
  
	phi.vec = c(1)
	for (ii in 1:length(i.data)){
		phi.vec = c(phi.vec, phi^ii)
	}
	exp.x = c(0)
	for (jj in 1:length(i.data)){
		cur.elements = (phi.vec[1:jj])*rev(i.data[1:jj])
		cur.expec = alpha * sum(cur.elements)
		exp.x = c(exp.x, cur.expec)
	}
	return(exp.x)
}

#Getting variance vector of xk , Var(xk)
V.x <- function(theta, n,sigma = sigma.fix) {
  phi = tanh(theta[2])
  fixed.var <- (sigma^2) / (1 - phi^2)
  var_vec = rep(fixed.var, (n + 1))
  
  return(var_vec)
}

#Tensors for states and parameters 

#Tensor for states, derived from expected Fisher information matrix 
tensorX   = function(beta = beta.fix, sigma = sigma.fix, theta, statevec, delta=1){
	mu = theta[1]
	phi = tanh(theta[2])
	alpha = theta[3]

	n = length(statevec)

	first  =  1 / (sigma^2)
	rest = (beta^2) * exp(mu + beta * statevec[2:(n-1)])*delta + (1+phi^2)/(sigma^2)
	last = (beta^2) * exp(mu + beta * statevec[n])*delta + 1/(sigma^2)
	diag.element = c(first,rest,last)

	supper = rep(-phi/(sigma^2), n-1)
    sub    = supper 

    tensor = triDiag(diag.element, supper, sub, nrow = n, ncol =n)

    return(tensor)
}

#Tensor for theta, derivded from expected Fisher information matrix 
tensorPar = function(beta = beta.fix, sigma = sigma.fix, theta, meanvec, varvec,i.data,delta = 1){
 	temp1 = 0
	mu = theta[1]
	phi = tanh(theta[2])
	alpha = theta[3]

	n = length(meanvec)

	#(1,1)
	temp1 = sum(exp(mu + beta * meanvec + (beta^2) * varvec *0.5) * delta)
	#for (ii in 1:length(meanvec)){
	#	cur.val = exp(mu + beta * meanvec[ii] + (beta^2) * varvec[ii]/2) * delta
	#	temp1 = temp1 + cur.val
	#}

    #(2,2) 
	temp2 = 2 * (phi^2) + n*(1-(phi^2)) + (1-(phi^2)) * sum(meanvec[1:(n-1)]^2) /(sigma^2)

	#(3,3)
	temp3 = sum(i.data^2) / (sigma^2)

	#(2,3) and (3,2)
	subsup = (1-(phi^2)) *sum(meanvec[1:(n-1)] * i.data) / (sigma^2)

	return(matrix(c(temp1, 0,0,0,temp2,subsup,0,subsup,temp3),nrow = 3, byrow = TRUE))
}


#Derivatives 

#Derivative of posterior w.r.t to X
Xder   = function(x.data,n.data,i.data,theta,beta = beta.fix, sigma = sigma.fix, delta = 1){
	mu = theta[1]
	phi = tanh(theta[2])
	alpha = theta[3]
	n = length(x.data)

	results = c()
	#k=0
	results[1] = (phi * x.data[2] - alpha * phi * i.data[1] - x.data[1])/(sigma^2)
	#k=1,...,k-1
	for (ii in 2:(n-1)){
		results[ii] =  (phi*x.data[ii+1] - (phi^2)*x.data[ii] - phi*alpha*i.data[ii]-x.data[ii] + phi*x.data[ii-1] + alpha*i.data[ii-1]) + beta * (n.data[ii-1]) - beta*exp(mu + beta * x.data[ii]) * delta
	}
	#k=K
	results[n] = (-x.data[n] + phi*x.data[n-1] + alpha * i.data[n-1])/(sigma^2) + beta*n.data[n-1] - beta*exp(mu + beta*x.data[n]) * delta 
	return(results)
}

#Derivative of posterior w.r.t parameter Theta 
Parder = function(x.data,n.data,i.data,theta,beta = beta.fix, sigma = sigma.fix, delta = 1){
	mu    = theta[1]
	phi   = tanh(theta[2])
	alpha = theta[3]

	n = length(x.data)

	Mu.der    = sum(n.data - exp(mu + beta * x.data[-1]) * delta)

	temp1 = sum(x.data[-n] * (x.data[-1] - phi * x.data[-n] - alpha * i.data ))
	Gamma.der = (1 - phi^2) * temp1 / (sigma^2) - phi + x.data[1]^2 * phi * (1-phi^2) / (sigma^2)

	temp2 = sum(i.data * (x.data[-1] - phi * x.data[-length(x.data)] - alpha * i.data))
	Alpha.der = temp2 / (sigma^2)

	return(c(Mu.der, Gamma.der, Alpha.der))
}

#Hamiltonians
#Hamiltonian for states 
Ham.x = function(x.data, n.data, i.data, theta, aux, Xten, Xten.inv){

	n = length(n.data)

	Ltemp     = -post.x(log = TRUE,x.data = x.data,n.data=n.data, i.data=i.data,theta= theta)
	logtemp   = 0.5 * log(((2*pi)^n) * det(Xten))
	potential = Ltemp + logtemp
	kinetic   = 0.5 * as.vector(t(aux)) %*% Xten.inv %*% as.vector(aux)

	return(potential + kinetic)
}

#Hamiltonian for parameter theta 
Ham.par = function(x.data, n.data, i.data, theta, aux, Parten,Parten.inv){
	

	Ltemp     = -post.theta(log = TRUE, x.data=x.data, n.data=n.data, i.data = i.data, theta = theta)
	logtemp   = 0.5 * log(((2*pi)^3) * det(Parten))
	potential = Ltemp + logtemp
	kinetic = 0.5 * as.vector(t(aux)) %*% Parten.inv %*% as.vector(aux)

	return(potential + kinetic)
}



#PDE's of the Hamiltonian for latent states 

#Partial derivatives of the hamiltonian of states w.r.t states   
Ham.Xderiv   = function(x.data, n.data, i.data,theta,aux,Xten,Xten.inv,beta=beta.fix,delta = 1){
	mu = theta[1]
	phi = tanh(theta[2])
	alpha = theta[3]

	#1st term
	postDeriv = Xder(x.data=x.data, n.data=n.data,i.data=i.data,theta=theta)

	#2nd
	matDeriv = c()
	matDeriv[1] = 0
	diag.temp = c(0,(beta^3) * exp(mu + beta * x.data[-1])*delta)

	for (ii in 2:length(x.data)){
		matDeriv[ii] = 0.5* diag(Xten.inv)[ii] * diag.temp[ii]
	}

	#3rd term 
	kinDeriv    = c()
  kinDeriv[1] = 0
  
	for(jj in 2:length(x.data)){

	  temp = sum(aux*Xten.inv[,jj]) * diag.temp[jj] * Xten.inv[jj,] 
    temp22 = temp %*% as.vector(aux)
		kinDeriv[jj] = 0.5 * temp22
	}

	return(as.vector(-postDeriv + matDeriv - kinDeriv))
}

#Partial derivative of the hamiltonian of theta w.r.t theta   
Ham.Parderiv = function(x.data,n.data,meanvec,varvec, i.data,theta,aux,Parten,Parten.inv,sigma=sigma.fix,delta = 1){
	mu = theta[1]
	phi = tanh(theta[2])
	alpha = theta[3]
	n = length(meanvec)

	#1st term 
	postDeriv = Parder(x.data = x.data, n.data = n.data, i.data=i.data, theta=theta)
  
	#Finding tensor derivative matrix 
	
	#For mu
	temp.mu = 0
	for (ii in 1:length(meanvec)){
	  cur.val = exp(mu + beta * meanvec[ii] + (beta^2) * varvec[ii]/2) * delta
	  temp.mu = temp.mu + cur.val
	}
	dGdmu = matrix(c(temp.mu,0,0,0,0,0,0,0,0),nrow = 3, byrow = TRUE)
	
	#For gamma
	temp1 = exp(mu + beta*meanvec + (beta^2) * varvec/2)
	temp2 = -beta * (sigma^2) * phi / ((1-phi^2)^2)
	der = c()
	der[1:2] = temp1[1:2] * temp2 * (1-phi^2)
	
	tempdir = c(1)
	for(jj in 3:(n-1)){
	  tempdir = c(tempdir,(jj-1)*phi^(jj-2))
	}
	
	dExdgamma = c(0,0)
	for(jj in 3:n){
	  dExdgamma[jj] = alpha*(sum(tempdir[1:(jj-2)] * rev(i.data[1:(jj-2)])))
	}
	
	for(jj in 3:n){
	  der[jj] = temp1[jj] * delta * (beta * dExdgamma[jj] + temp2)
	}

	
	temp.11 = (1-phi^2) * sum(der) 
	temp.22 = (4*phi-2*n*phi-4*phi*sum(meanvec[1:(n-1)] * dExdgamma[1:(n-1)]) /(sigma^2)) * (1-phi^2)
	temp.23 = temp.32= -2*phi *sum(dExdgamma[1:(n-1)] * i.data) * (1-phi^2) / (sigma^2)
	
	dGdgamma = matrix(c(temp.11,0,0,0,temp.22,temp.23,0,temp.32,0),nrow=3,byrow=TRUE)

	#For alpha
	phi.vec = c(1)
	for (ii in 1:length(i.data)){
	  phi.vec = c(phi.vec, phi^ii)
	}

	dExdalpha = c(0)
	for (jj in 1:length(i.data)){
	  cur.elements = (phi.vec[1:jj])*rev(i.data[1:jj])
	  cur.expec = sum(cur.elements)
  	  dExdalpha = c(dExdalpha, cur.expec)
	}
  
	tempalpha=c()
	for(kk in 1:n){
	  tempalpha[kk] = temp1[kk] * delta * (beta *dExdalpha[kk])
	}
	temp.11.alpha = sum(tempalpha)
	temp.22.alpha = 2*(1-phi^2)*sum(meanvec[1:(n-1)] * dExdalpha[1:(n-1)]) / (sigma^2)
	temp.23.alpha = temp.32.alpha = (1-phi^2)*sum(meanvec[1:(n-1)] * dExdalpha[1:(n-1)]) / (sigma^2)

	dGdalpha = matrix(c(temp.11.alpha,0,0,0,temp.22.alpha,temp.23.alpha,0,temp.32.alpha,0),nrow=3,byrow= TRUE)
	
	#2nd term
  deriv.mu    = 0.5 * sum(diag((Parten.inv %*% dGdmu)))
	deriv.gamma = 0.5 * sum(diag((Parten.inv %*% dGdgamma)))
  deriv.alpha = 0.5 * sum(diag((Parten.inv %*% dGdalpha)))


	#3rd term 
	mu.temp3    = 0.5 * as.vector(aux) %*% Parten.inv %*% dGdmu %*% Parten.inv %*% as.vector(aux)
	gamma.temp3 = 0.5 * as.vector(aux) %*% Parten.inv %*% dGdgamma %*% Parten.inv %*% as.vector(aux)
  alpha.temp3 = 0.5 * as.vector(aux) %*% Parten.inv %*% dGdalpha %*% Parten.inv %*% as.vector(aux)

	matDeriv = c(deriv.mu, deriv.gamma, deriv.alpha)
	kinDeriv = c(mu.temp3, gamma.temp3, alpha.temp3)

	return(as.vector(-postDeriv + matDeriv - kinDeriv))
}



#Partial derivative of Hamiltonians with respect to auxiliary variables
Ham.Aux.deriv = function(aux,inv){
	return(inv %*% as.vector(aux))
}


#Use leapfrog method to simulate a trajectory of Hamiltonian system
gen.X = function(x.data, n.data,i.data,theta, aux, Xten, Xten.inv,step.size = 0.02, steps = 25){
	while(steps>0){
	  #print(paste("Current running",26-steps,"in Leapfrog integrator"))
		(aux.half = aux - 0.5 * step.size * Ham.Xderiv(x.data = x.data, n.data = n.data, i.data = i.data, theta=theta,aux = aux,Xten = Xten, Xten.inv = Xten.inv))
		(x.data   = x.data +  step.size * Ham.Aux.deriv(aux = aux.half, inv = Xten.inv))
		(aux      =  aux.half - 0.5*step.size * Ham.Xderiv(x.data =x.data,n.data=n.data,i.data=i.data, theta=theta,aux=aux.half, Xten = Xten, Xten.inv=Xten.inv))
		steps = steps - 1
	}

	return(list(Xcandidate = x.data, aux = aux))
}

#Use leapfrog method to simulate a trajectory of Hamiltonian system
gen.Par = function(x.data,n.data,i.data,meanvec,varvec,theta,aux,Parten,Parten.inv,step.size = 0.05, steps = 15){
	while(steps>0){
	  #print(paste("Current running",26-steps,"in Leapfrog integrator"))
		aux.half = aux - 0.5 * step.size * Ham.Parderiv(x.data=x.data,n.data=n.data,i.data=i.data,meanvec = meanvec,varvec = varvec,theta=theta,aux = aux, Parten = Parten,Parten.inv = Parten.inv)
		theta = theta + step.size * Ham.Aux.deriv(aux = aux.half,inv = Parten.inv)
		aux   = aux.half - 0.5*step.size * Ham.Parderiv(x.data=x.data,n.data=n.data,i.data=i.data,meanvec = meanvec,varvec = varvec,theta=theta,aux = aux.half, Parten = Parten.inv,Parten.inv = Parten.inv)

		steps = steps - 1
	}

	return(list(Parcandidate = theta, aux = aux))
}


