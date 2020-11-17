cat("Reading BM.R\n")


# global variables:
nSteps <- 24
h <- 0.25                            # stepsize
d <- seq(-1,-4,length.out=nSteps)    # death thresholds


# H(j) from docs/BrownianMotionProblem.pdf, eq(2)
#
fcn_H <- function(j){
   
   t_j <- j*h
   sqrt_tj <- sqrt(t_j)
   sqrt_h <- sqrt(h)
   integrand <- function(y)
      pnorm((y-d[j+1])/sqrt_h)*dnorm(y,mean=0,sd=sqrt_tj)
   integrate(integrand,
      lower = d[j], upper = Inf, rel.tol = 1e6*.Machine$double.eps
   )$value / pnorm(-d[j]/sqrt_tj)
}


# H(j) from docs/BrownianMotionProblem.pdf, p2, eqa(4),(5).
# This time implemented using multivariate normal cumulative
# distribution function from package mvtnorm.
#
fcn_H1 <- function(j){
   
   # Cov(Y(t_j)>=d_j,Y(t_{j+1}))
   tj <- j*h
   C <- rbind(c(tj,tj),c(tj,tj+h))
   # P(Y(t_j)>=d_j and Y(t_{j+1})>=d_{j+1})
   mu <- rep(0,2)
   lb <- c(d[j],d[j+1])
   ub <- c(Inf,Inf)
   num <- mvtnorm::pmvnorm(lower=lb, upper=ub,mean=mu,sigma=C)
   # P(Y(t_j)>=d_j)
   denom <- pnorm(-d[j]/sqrt(tj)) 
   num/denom
}


# Survival probabilities p_j, 
# see docs/BrownianMotionProblem.pdf, eqs(2),(3)
#
survivalProbs <- function(){
   
   p <- rep(0,nSteps)
   p[1] <- pnorm(-d[1]/sqrt(h))
   for(j in 1:(nSteps-1)) p[j+1] <- p[j]*fcn_H(j)
   p
}
# Death probabilities q_j, 
# see docs/BrownianMotionProblem.pdf, eqs(2),(3)
#
deathProbs <- function(){
   
   # probability of survival to time t=0 is one:
   p <- c(1,survivalProbs())
   vapply(1:nSteps,FUN=function(j) p[[j]]-p[[j+1]],FUN.VALUE=0)
}


# Path of a standard Brownian motion (starting at zero) evaluated at time points
# t_j = j*h, j=0,1,..,nSteps.
#
nextBrownianPath <- function()
   c(0,cumsum(sqrt(h)*rnorm(nSteps)))


# pValue of deviation (2-sided) of sample mean from theoretical mean.
# Distribution of sample mean approximated by normal distribution, i.e.:
#
#      P[abs(X-mu)>=abs(realizedSampleMean-mu)]
#
# where X=(X_1+X_2+...+X_n)/n is the sample mean.
# See docs/BrownianMotionProblem.pdf, p2, eq(8).
#
# @param mu theoretical mean (population mean)
# @param sigma theoretical sd (population sd)
#
pValueOfSampleMeanDeviation <- function(mu,sigma,sampleMean,sampleSize)
   2*pnorm(-sqrt(sampleSize)*abs(sampleMean-mu)/sigma)



# Simulation of Brownian paths, counting deaths at each time
# t_j = j*h.
# Results gathered in data.frame and written to file results/DeathStats.txt
#
runSimulation <- function(nPaths){
   
   cat("\n\nRunning simulation...")
   q <- deathProbs()      # probability of death at each time t_j
   Edc <- nPaths*q        # expected death counts at each time t_j
   dc <-  rep(0,nSteps)   # realized death counts at each time t_j
   
   
   for(i in 1:nPaths){
      
      Y <- nextBrownianPath() # starts with Y(0)=0, length 1+nSteps
      # step along path, check if death occurs
      j <- 1                  # time step
      dead <- FALSE
      while(j<=nSteps & !dead){
         
         if(Y[j+1]<d[j]){ dc[[j]] <- dc[[j]]+1; dead <- TRUE }     # death at time t_j recorded
         j <- j+1
      }
   }
   # p-values of 2-sided deviation of realized from expected death counts
   # note that counts are sums of indicator variables (0/1 valued), thus empirical probability
   # of death is a sample mean
   pVals <- vapply(
      1:nSteps,
      FUN=function(j){ 
         sigma <- sqrt(q[j]*(1-q[j]))    # sd of indicator variable P(X=1)=p[j]=1-P(X=0)
         pValueOfSampleMeanDeviation(q[j],sigma,dc[[j]]/nPaths,nPaths)
      },
      FUN.VALUE=0
   )
   result <- data.frame(
      Step = 1:nSteps,
      DeathCountsActual = dc,
      DeathCountsExcpected = round(Edc),
      DeviationPcnt = round(100*(dc-Edc)/Edc,2),
      pValue = round(pVals,4)
   )
   outFile <- "results/DeathStats.txt"
   capture.output(print(result),file=outFile)
   cat("Finished, results in",outFile,"\n")
   (result)
}
