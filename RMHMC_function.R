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
prior.x = function(log = TRUE, x.data){
  priors = sum(dnorm(x.data,mean = 0 , sd = 10, log = TRUE))
  if (log) {
  	return(priors)
  }else{
  	return(exp(priors))
  }
}
#flat priors of theta 
prior.par = function(log = TRUE, theta){
	priors = sum(dnorm(theta, mean = 0 , sd = 10, log = TRUE))
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
	
	for (ii in 1:length(meanvec)){
		cur.val = exp(mu + beta * meanvec[ii] + (beta^2) * varvec[ii]/2) * delta
		temp1 = temp1 + cur.val
	}
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
	
	results = c()
	results[1] = (phi * x.data[2] - alpha * phi * i.data[1] - x.data[1])/(sigma^2)
	for (ii in 2:length(x.data)){
		results[ii] =  -(x.data[ii] - phi * x.data[ii-1] - alpha*i.data[ii-1])/(sigma^2) + beta * (n.data[ii-1] - exp(mu + beta * x.data[ii-1]) * delta)
	}
	return(results)
}

#Derivative of posterior w.r.t parameter Theta 
Parder = function(x.data,n.data,i.data,theta,beta = beta.fix, sigma = sigma.fix, delta = 1){
	mu    = theta[1]
	phi   = tanh(theta[2])
	alpha = theta[3]

	Mu.der    = sum(n.data - exp(mu + beta * x.data[-1]) * delta)

	temp1 = sum(i.data * (x.data[-1] - phi * x.data[-length(x.data)] - alpha * i.data))
	Alpha.der = temp1 / (sigma^2)

	temp2 = sum(x.data[-length(x.data)] * (x.data[-1] - phi * x.data[-length(x.data)] - alpha * i.data ))
	Gamma.der = (1 - phi^2) * temp2 / (sigma^2) - phi + x.data[1]^2 * phi * (1-phi^2) / (sigma^2)

	return(c(Mu.der, Alpha.der, Gamma.der))
}

#Hamiltonians
#Hamiltonian for states 
Ham.x = function(x.data, n.data, i.data, theta, aux, Xten, Xten.inv){

	D = length(n.data)

	Ltemp     = -post.x(log = TRUE,x.data = x.data,n.data=n.data, i.data=i.data,theta= theta)
	logtemp   = 0.5 * (D*log(2*pi) + determinant(Xten,logarithm=TRUE)$modulus)
	potential = Ltemp + logtemp
	kinetic   = 0.5 * as.vector(t(aux)) %*% Xten.inv %*% as.vector(aux)

	return(potential + kinetic)
}

#Hamiltonian for parameter theta 
Ham.par = function(x.data, n.data, i.data, theta, aux, Parten,Parten.inv){
	

	Ltemp     = -post.theta(log = TRUE, x.data=x.data, n.data=n.data, i.data = i.data, theta = theta)
	logtemp   = 0.5 * (3 * log(2*pi) + determinant(Parten,logarithm=TRUE)$modulus)
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

	diag.temp = diag(Xten.inv)
	for (ii in 2:length(x.data)){
		matDeriv[ii] = 0.5* diag.temp[ii] * (beta^3) * exp(mu + beta * x.data[ii]) * delta 
	}

	#3rd term 
	kinDeriv    = c()
	kinDeriv[1] = 0

	for(jj in 2:length(x.data)){
		kinDeriv[jj] = sum(aux* Xten.inv[,jj]) *  (beta^3) * exp(mu + beta * x.data[jj]) * delta * diag.temp[jj] * aux[jj]
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

	#2nd term
	#for mu 
	temp.mu = 0
	#(1,1)
	for (ii in 1:length(meanvec)){
		cur.val = exp(mu + beta * meanvec[ii] + (beta^2) * varvec[ii]/2) * delta
		temp.mu = temp.mu + cur.val
	}
	deriv.mu = 0.5 * Parten.inv[1,1] * temp.mu 

	#for gamma
	temp.22 = 4*phi-2*n*phi-2*(1-phi)*sum(meanvec[1:(n-1)]^2) * (1-phi^2)/(sigma^2)
	temp.23 = -2*phi *sum(meanvec[1:(n-1)] * i.data) * (1-phi^2) / (sigma^2)
	mat.temp = matrix(c(0,0,0,0,temp.22,temp.23,0,temp.23,0),nrow=3)
	deriv.gamma = 0.5 * sum(diag((Parten.inv %*% mat.temp)))

	#for alpha 
	#dH_da = 0

	#3rd term 
	#for mu 
	mu.temp3 = 0.5* sum(aux*Parten.inv[,1]) * temp.mu * Parten.inv[1,1] * aux[1]
	gamma.temp3 = 0.5 * as.vector(aux) %*% Parten.inv %*% mat.temp %*% Parten.inv %*% as.vector(aux)


	matDeriv = c(deriv.mu, deriv.gamma, 0)
	something = c(mu.temp3, gamma.temp3, 0)

	return(as.vector(-postDeriv + matDeriv - something))
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
gen.Par = function(x.data,n.data,i.data,meanvec,varvec,theta,aux,Parten,Parten.inv,step.size = 0.02, steps = 20){
	while(steps>0){
	  #print(paste("Current running",26-steps,"in Leapfrog integrator"))
		aux.half = aux - 0.5 * step.size * Ham.Parderiv(x.data=x.data,n.data=n.data,i.data=i.data,meanvec = meanvec,varvec = varvec,theta=theta,aux = aux, Parten = Parten,Parten.inv = Parten.inv)
		theta = theta + step.size * Ham.Aux.deriv(aux = aux.half,inv = Parten.inv)
		aux   = aux.half - 0.5*step.size * Ham.Parderiv(x.data=x.data,n.data=n.data,i.data=i.data,meanvec = meanvec,varvec = varvec,theta=theta,aux = aux.half, Parten = Parten.inv,Parten.inv = Parten.inv)

		steps = steps - 1
	}

	return(list(Parcandidate = theta, aux = aux))
}


