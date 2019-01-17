#### (invisible function)
#### Approximation of the negative log-gamma distribution by Gaussian mixtures
#### Compute mixture parameters, given degrees of freedom (n)
#### with R(v) components, weights w(v), means m(v) and variances v(v) 

## This file is based on the bayesf MATLAB package found on S. Fruehwirth-Schnatters'
## website 
## <http://statmath.wu.ac.at/~fruehwirth/monographie/>.

#### last version: 2019/01/07
#### This function is not meant to be called directly by the user. 

#------ individual mixtures -------------------------------------------------- -
# R(v)= 10 components for 1 <= v <= 4         ... range = 1 
# R(v)= 9 components for 5 <= v <= 19         ... range = 2  
#----- parameterized mixtures ------------------------------------------------ -
# R(v)= 4 components for 20 <= v <= 49        ... range = 3
# R(v)= 3 componnets for 50 <= v <= 439       ... range = 4
# R(v)= 2 components for 440 <= v <= 1599     ... range = 5
# R(v)= 2 components for 1600 <= v <= 10000   ... range = 6
# R(v)= 2 components for 10000 <= v <= 30000  ... range = 7
#----------------------------------------------------------------------------- -

compute_mixture <- function(n, mcomp){
  
  # get individual or parameterized mixture components from mixcomp_poisson()
  if(n <= 30000){
    m <- mcomp$m[n, ]
    v <- mcomp$v[n, ]
    w <- mcomp$w[n, ]
  } else {
    m <- 0
    v <- 1
    w <- 1
    
    m <- sqrt(trigamma(n))*m + (-1)*digamma(n)
    v <- trigamma(n)*v
  }
  
  nc <- length(m)
  out <- list(m = m, v = v, w = w)
  return(out)
}

