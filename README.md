# BM
## Brownian motion hitting a lower barrier

Fix a step size h > 0 and values 
   
     0 > d_0 > d_1 > ... > d_n
    
(lower barrier values). A standard Brownian motion Y(t) is sampled at equidistant times t_j = j*h    
yielding the process Y_j=Y(t_j) which is a simple random walk with IID increments distributed as
sqrt(h)*N(0,1).

As such Y_j is a Markov process. We say that the process Y_j dies at time t_j if

   Y_k >= d_k, for k<j, but Y_j < d_j
   
i.e. the process Y_j dies as soon as it falls below the lower barrier d.
Now we compute the probability p_j that Y_j is still alive at time t_j by using the following intuitive
application of the Markov property:
 
     P( Y_{j+1}>=d_{j+1} | Y_j>=d_j, Y_{j-1}>=d_{j-1},...,Y_1>=d_1 ) = P( Y_{j+1}>=d_{j+1} | Y_j>=d_j ).   (1)
     
Sloppy intuitive formulations of the Markov property suggest as much but indeed the Markov property does not
imply this and in our (extremely simple) case it is in fact not true.
 
This leads to incorrect values for the survival probabilities p_j.
We also compute the survival pobabilities p_j correctly (directly from the multivariate normal distribution function
obtained from the R-package mvtnorm) and compute the deviation from the incorrect values.

The correct theoretical values are then checked via simulation.

This project came about from an effort to debug of a more complicated project, where the Markov property 
was incorrectly applied as in (1).
Now it only serves as a graphic demonstration of the failure of (1) for the simplest Markov processes.
