## multipar EVPI via GP
## Markov model case study
## June 2012

library(MASS) # for mvrnorm
library(gtools) # for rdirichlet

theta0 <- c(1000.00, 0.10, 5.20, 400.00, 0.70, 0.30, 3.00, 0.25, 
	-0.10, 0.50, 1500.00, 0.08, 6.10, 0.80, 0.30, 3.00, 0.20, -0.10, 0.50)

stdvector<-c(1, .02, 1, 200, 0.1, 0.1, .5,0.1, .02, .2, 1, .02, 1, .1,.05,  1, .05, .02, .2)
n.inputs<-19

rho <- 0.6
cormatrix<-matrix(0, n.inputs, n.inputs)
diag(cormatrix) <- 1

cormatrix[5,7] <- cormatrix[7,5] <- rho
cormatrix[5,14] <- cormatrix[14,5] <- rho
cormatrix[5,16] <- cormatrix[16,5] <- rho
cormatrix[7,14] <- cormatrix[14,7] <- rho
cormatrix[7,16] <- cormatrix[16,7] <- rho
cormatrix[14,16] <- cormatrix[16,14] <- rho

sigma0 <- diag(stdvector) %*% cormatrix %*% diag(stdvector)

lambda <- 10000

model.func <- function(params, lambda, T) {   
	p <- params[1:19]
	M0 <- matrix(params[20:28], nrow=3, byrow=T)
	M1 <- matrix(params[29:37], nrow=3, byrow=T)

    response <- rep(0, len = T)

    S0 <- matrix(c(p[5], 1 - p[5], 0), nrow = 1)
    U0 <- matrix(c(p[6], 0, 0), nrow = 3)
    response0 <- S0 %*% U0
    SM00 <- S0 %*% M0

    S1 <- matrix(c(p[14], 1 - p[14], 0), nrow = 1)
    U1 <- matrix(c(p[15], 0, 0), nrow = 3)
    response1 <- S1 %*% U1
    SM01 <- S1 %*% M1

    for(i in 1:T) {
        response0[i] <- SM00 %*% U0
        SM00 <- SM00 %*% M0
        response1[i] <- SM01 %*% U1
        SM01 <- SM01 %*% M1
    }
    nb0<-lambda * (sum(response0,response) + p[8] * p[9] * p[10]) - (p[1] + p[2] * p[3] * p[4])
    nb1<-lambda * (sum(response1,response) + p[17] * p[18] * p[19]) - (p[11] + p[12] * p[13] * p[4])
	return(c(nb0, nb1))
}


gen.inputs<- function(ndraw,theta0, sigma0, lambda=lambda, T=T) {
	theta <- MASS::mvrnorm(ndraw, mu = theta0, Sigma = sigma0)

	responding0 <- gtools::rdirichlet(ndraw,c(7, 4, 2))
	notresponding0 <- gtools::rdirichlet(ndraw,c(1, 10, 2))

	responding1 <- gtools::rdirichlet(ndraw,c(7, 4, 2))
	notresponding1 <- gtools::rdirichlet(ndraw,c(1, 10, 2))

	died <- matrix(rep(c(0, 0, 1), ndraw), nrow=ndraw, byrow=T)

	M0 <- cbind(responding0, notresponding0, died)
	M1 <- cbind(responding1, notresponding1, died)

	cbind(theta, M0, M1)  
}

s <- 10000
T <- 20

set.seed(1)
input.mat <- gen.inputs(ndraw=s, theta0=theta0, sigma0=sigma0, lambda=lambda, T=T)

nb <- t(apply(input.mat, 1, model.func, lambda, T))

colMeans(nb)
max(colMeans(nb))
mean(apply(nb, 1, max))
mean(apply(nb, 1, max)) - max(colMeans(nb))

# parallelise

# first allow model function to be applied over a list
input.mat.df <- data.frame(t(input.mat)) # define inputs as data.frame (i.e. list) object
nb.p <- t(simplify2array(lapply(input.mat.df, model.func, lambda, T)))

colMeans(nb.p)
max(colMeans(nb.p))
mean(apply(nb.p, 1, max))
mean(apply(nb.p, 1, max)) - max(colMeans(nb.p))


library(parallel)
mc.cores <- 4 # assume 4 cores

input.mat.df <- data.frame(t(input.mat))
nb.p <- t(simplify2array(mclapply(input.mat.df, FUN = model.func, mc.cores = mc.cores, lambda, T)))

colMeans(nb.p)
max(colMeans(nb.p))
mean(apply(nb.p, 1, max))
mean(apply(nb.p, 1, max)) - max(colMeans(nb.p))


