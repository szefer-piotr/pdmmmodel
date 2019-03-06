# $Id: nichespace.R,v 1.3.1.1 2010/05/04 17:45:58 axel Exp axel $

# (c) 2008 by Axel G. Rossberg and Leo Nagelkerke

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# SUMMARY:

# This is experimental code to estimate the function a(t,s) that
# determines the interaction strength between two species (or
# individuals) based on vectors of traits of "consumer" (t) and
# "resource" (s).

# niche.quadratic.form(...): estimate quadratic form from empirical data
# niche.mock.data(...)     : generate mock data for testing purposes
# niche.a(...)             : compute interaction strength a

################################################################
################################################################

# This is an implementation of a scalar (dot) product between vectors.
# There must be a different definition of a scalar product in R!
# Could anybody help?

niche.dot <- function (x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  return(as.double(t(x) %*% y))
}

################################################################
################################################################

# Compute interaction strength a from vectors of resource traits
# traits.r, consumer traits traits.c, and a precomputed quadratic form
# qf, generated, for example, by niche.quadratic.form(...)

niche.a <- function(traits.r,traits.c,qf){
  return(exp(niche.loga(traits.r,traits.c,qf)))
}

niche.loga <- function(traits.r,traits.c,qf){
  # Compute pair traits:
  v <- c(traits.r,traits.c)
  # Evaluate quadratic form:
  return(qf$a + niche.dot(qf$b, v) + 0.5*niche.dot(v, qf$C %*% v) )
}

################################################################
################################################################

# This function generates random mock data for testing purposes.

## Arguments: ..
##  
## n.r : number of resource species (default 20)
##
## n.c : number of consumer species (default 30)
##
## ntraits.r : number of phenotypic traits determined for resources (def. 2)
##
## ntraits.c : number of phenotypic traits determined for consumer (def. 3)
##
## b : the vector b of the quadratic form or a scalar taken as the
##     value for all components of b (default 0)
##
## C : the matrix C of the quadratic form.  The default is chosen such
##     that the quadratic part of the quadratic form evaluates to
##     (-5/2)*(traits.r(1)-traits.c(1))^2, i.e., such that only the
##     leading components of the trait vectors are trophic traits.

## Output: ..
##
## a list with named elements traits.r, traits.c, interactions, which
## can serve as the corresponding input data for niche.quadratic.form
  
niche.mock.data <- function (n.r=20,
                             n.c=30,
                             ntraits.r=2,
                             ntraits.c=3,
                             b=0,
                             C=NA
                             ){

  # Generate b vector depending on input argument:
  b <- as.vector(array(b,dim=ntraits.c+ntraits.r))

  # Generate C matrix depending on input argument:
  if(is.na(C)){
    C=diag(0,ntraits.c+ntraits.r);
    C[          1,          1]=-5
    C[ntraits.r+1,ntraits.r+1]=-5
    C[1          ,ntraits.r+1]=5
    C[ntraits.r+1,          1]=5
  }

  # Generate random traits for resources and consumers:
  t=array(rnorm(n.r*ntraits.r),dim=c(n.r,ntraits.r))
  s=array(rnorm(n.c*ntraits.c),dim=c(n.c,ntraits.c))

  # Compute interaction strengths between resources and consumers
  a <- array(dim=c(n.r,n.c))
  for(i in 1:n.r){
    for(j in 1:n.c){
      v <- c(t[i,],s[j,])
      a[i,j] <- exp(niche.dot(b,v) + 0.5 * niche.dot(v,C %*% v))
    }
  }

  return(list( traits.r=t, traits.c=s, interactions=a ))
}


################################################################
################################################################

# This function estimates the quadratic form from traits and
# interaction strengths of specific species, taking optional
# population abundances and statistical weights into account.

## Arguments: ..
## 
## traits.r : matrix containing in each row the traits of one resource
## 
## traits.c : matrix containing in each row the traits of one consumer
## 
## interactions : matrix containing the strengths of observed
##                interactions between resources (rows) and consumers
##                (columns)
## 
## weights.r : vector of weights for resource species (default: all equal)
## 
## weights.c : vector of weights for consumer species (default: all equal)
## 
## rho.r : vector of abundances of resource species (default: all equal)
## 
## rho.c : vector of abundances of resource species (default: all equal)

## Output: ...
## 
## a list with named elements b, C, and a. Here b and C are the
## ESTIMATED coefficients of the quadratic form, and a is a function
## that computes the ESTIMATED interaction strength (up to a factor)
## from two trait vectors as arguments.

niche.quadratic.form <- function (traits.r,traits.c,interactions,
                                  weights.r=1,
                                  weights.c=1,
                                  rho.r=1,
                                  rho.c=1
                                  ){
  # Introduce abbreviations t and s in line with theory:
  t <- as.matrix(traits.r)
  s <- as.matrix(traits.c)

  # Some more abbreviations:
  n.r <- dim(t)[1]
  n.c <- dim(s)[1]
  ntraits.r <- dim(t)[2]
  ntraits.c <- dim(s)[2]
  
  # Make sure weights and abundances are of the right form, and
  # introduce abbreviations g and h in line with theory:
  g <- array(as.vector(weights.r),n.r)
  h <- array(as.vector(weights.c),n.c)
  rho.r <- array(as.vector(rho.r),n.r)
  rho.c <- array(as.vector(rho.c),n.c)
  
  # Normal approximation for resource trait distribution:
  C.u <- cov.wt(t, wt=as.vector(rho.r*g), method="ML")
  # (Note: Since we are using weights, "ML" is the better choice!)
  mu.u <- C.u$center
  C.u <- C.u$cov

  # Normal approximation for consumer trait distribution:
  C.v <- cov.wt(s, wt=as.vector(rho.c*h), method="ML")
  mu.v <- C.v$center
  C.v <- C.v$cov
  
  # Compute an array which contains in each row the traits of a
  # consumer-resource pair:

  pair.traits <- array(dim=c(n.r*n.c,ntraits.r+ntraits.c))
  colnames(pair.traits) <- c(colnames(t),colnames(s))

##   print(dim(pair.traits))
##   print(dim(t))
##   print(dim(s))
##   print(n.r)
##   print(n.c)
  
  for(i in 1:n.r){
    for(j in 1:n.c){
      pair.traits[(j-1)*n.r+i,] <- c(t[i,],s[j,])
    }
  }
  pair.weights <- as.vector(interactions * outer(g,h))
  
  # Normal approximation for pair-traits observed in interactions:
  C.p <- cov.wt(pair.traits, wt=pair.weights, method="ML")
  mu.p <- C.p$center
  C.p <- C.p$cov

  # Invert covariances
  C.u.inv <- solve(C.u)
  C.v.inv <- solve(C.v)
  C.p.inv <- solve(C.p)

  # Compute b
  b <- C.p.inv %*% mu.p - rbind(C.u.inv %*% mu.u,C.v.inv %*% mu.v)

  # Compute C
  C1 <- -C.p.inv
  C2 <- cbind(rbind(C.u.inv,array(0,dim=c(ntraits.c,ntraits.r))) ,
              rbind(array(0,dim=c(ntraits.r,ntraits.c)),C.v.inv) )
  C <- C1 + C2

  # Define interaction strength function
  afunc <- function(traits.r,traits.c){
    v <- c(traits.r,traits.c)
    return(exp( niche.dot(b, v) + 0.5*niche.dot(v, C %*% v) ))
  }

  unscale.interactions <- NULL
  
  for(c in 1:n.c){
    for(r in 1:n.r){
      unscale.interactions <- c(unscale.interactions,
                                afunc(t[r,],s[c,])*
                                rho.r[r]*rho.c[c]*g[r]*h[c]
                                )
    }
  }

  
  exp.a <- mean(pair.weights)/mean(unscale.interactions)
  a <- log(exp.a)
  
  mu.prime <- -solve(C,b)
  a.prime <- a+niche.dot(mu.prime,b)/2
  
  return( list(a=a, b=b, C=C, a.prime=a.prime, mu.prime=mu.prime ))
}

#### Examples:
##
# mock <- niche.mock.data() # compute random mock data
# 
# QF <- niche.quadratic.form(mock$traits.r,mock$traits.c,mock$interactions)
# 
# eigen(QF$C) # There is one dominant eigenvalue corresponding to a
            # single pair of trophic traits.  The corresponding
            # eigenvector (last column) shows which traits are involved.

# QF$a(c(1,2),c(0,2,1))  # compute an interaction strength
# 
# niche.a(c(1,2),c(0,2,1),QF) # the same
