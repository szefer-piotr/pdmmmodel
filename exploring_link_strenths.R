rm(list = ls())
source("nagelkerke_rossberg_suppl/nagelkerke_rossberg_Rcode/nichespace.R")

# Set some parameters for the random community ------
n.r = 25
n.c = 36
ntraits.r=4
ntraits.c=6

# Random normal trait data
mock <- niche.mock.data(n.r = 25,
                        n.c = 36,
                        ntraits.r=4,
                        ntraits.c=6,
                        b=0,
                        C=NA)

QF <- niche.quadratic.form(mock$traits.r,mock$traits.c,mock$interactions)


# Extracting vectors and matrices ----------------------
# What is x? x is any combination of resource and consumer 
# traits concatenated together
# would that be all cobinations of resorce consumer?
# Yes, because it would have to fill the whole interaction matrix

# One random species
x <- c(mock$traits.r[sample(1:n.r, 1),], mock$traits.c[sample(1:n.c,1),])

b <- QF$b
C <- QF$C
a <- QF$a
u <- eigen(QF$C)$vectors # There is one dominant eigenvalue corresponding to a
# single pair of trophic traits.  The corresponding
# eigenvector (last column) shows which traits are involved.

lk <- eigen(QF$C)$values

# Traits vuln + foraging
Dprime = ntraits.c + ntraits.r
v <-  x[1:ntraits.r]

# Testing the model
# Logarithm of the link strength function
a + (t(b) %*% x) + 0.5 * (t(x) %*% (C %*% x))
niche.loga(v,x[(ntraits.r+1):(ntraits.r+ntraits.c)],QF)

# Link strength
exp(a + niche.dot(b, x) + 0.5 * niche.dot(x, C %*% x))
niche.a(x[1:ntraits.r],x[(ntraits.r+1):(ntraits.r+ntraits.c)],QF)


# Equation 8.11 Reconstruct the C matrix  ----------
C_rec <- matrix(0, nrow = Dprime, ncol = Dprime)

for(k in 1:Dprime){
  C_k <- (as.matrix(u[,k])*lk[k]) %*% t(as.matrix(u[,k]))
  C_rec <- C_rec + C_k
}

# This is true up to some rounding
round(C,4) == round(C_rec, 4)


# Equation 8.12 ------------
# for species 1

species = 1
x <- c(mock$traits.r[species,], 
       mock$traits.c[species,])

v = x[1:ntraits.r]

f = x[(ntraits.r+1):(ntraits.r+ntraits.c)]

# What about x = sum(u_k (u_k^Tx))
x_rec <- vector("numeric", length = Dprime)

for (k in 1:Dprime){
  x_k <- u[,k] %*% (t(u[,k]) %*% x)
  x_rec <- x_rec + x_k
}

round(x_rec,4) == round(x, 4)

# Equation 8.13 after decomposition ---------------------------

buux <- 0
for (k in 1:Dprime){
  buux <- buux + (t(b)%*%u[,k] * t(u[,k])%*%x)
}

lkux2 <- 0
for (k in 1:Dprime){
  lkux2 <- lkux2 + lk[k] * (t(u[,k])%*%x)^2
}

# log link
ll <- a + buux + 0.5 * lkux2

exp(ll)

# compare with a value predicted by function niche.a
niche.a(v,f,QF)


## Reduction of the number of trophic traits ------------------

# If  the $u^T_kX$ are approximately normally distributed
# as they likely will be because they are given by the sums 
# of independent random numbers. 
# $$ \frac{1}{4} \sum_{k=D+1}^{D^{'}} var[\lambda_k(\mathbf{u}^T_k\mathbf{x})^2 ] = \frac{1}{2} \sum_{k=D+1}^{D^{'}}\labda^2_k $$
# Choosing the truncation D large enough that this erroe remains within 
# acceptable bounds:

## Equation 8.15 -----------------
D <- 10

lk_inv = lk[order(lk)]
u_inv = u[, order(lk)]

buux_Dprime <- 0
# for (k in 1:Dprime){
#   buux_Dprime <- buux_Dprime + (t(b)%*%u_inv[,k] * t(u_inv[,k])%*%x)
# }

for (k in 1:Dprime){
  buux_Dprime <- buux_Dprime + (t(b)%*%u[,k] * t(u[,k])%*%x)
}

lkux2_D <- 0

for (k in 1:D){
  lkux2_D <- lkux2_D + lk[k] * (t(u[,k])%*%x)^2
}

# for (k in 1:D){
#   lkux2_D <- lkux2_D + lk_inv[k] * (t(u_inv[,k])%*%x)^2
# }

lk_2 <- sum((lk[(D+1):Dprime])/2)
lk_2 <- 0
# log link

llreal <- a + buux_Dprime + 0.5 * lkux2_D + lk_2

exp(llreal)


# Why it does not give a better approximation??

# Are eigenvalues not ordered????

llappr <- a + buux + 0.5 * lkux2_D + lk_2

exp(llappr)

## Why it doesnt match the real value?
niche.a(v,f,QF)



## Equation 8.16 ----------------------------------------------

t(u[,1]) %*% as.matrix(x)

v <- mock$traits.r
f <- mock$traits.c

g1v <- t(u[1:ntraits.r, 1]) %*% v[1, ]
h1f <- t(u[(ntraits.r+1):(ntraits.r + ntraits.c), 1]) %*% f[1, ] 

Vk <- g1v + h1f
Fk <- -h1f

# Baseline trophic traits
Vbase <- 0
for(k in 1:Dprime){
  Vbase <- Vbase + t(u[1:ntraits.r, k]) %*% v[k, ]
}

Fbase <- 0
for(k in 1:Dprime){
  Fbase <- Fbase + t(u[(ntraits.r+1):(ntraits.r + ntraits.c), k]) %*% f[k, ] 
}

# Equation 8.19 -----------------------------------------------

lkVkFk <- 0
for ()

logA = a + (Vbase + Fbase + 0.5 * lkVkFk)

















# Components of v and f have to be uncorrelated and have variance 1
# and mean 0, effect of their concatenation should therefore have 
# the same properties.

vec <- c()
for (k in 1:Dprime){
  vec <- c(vec, (t(u[,k])%*%x))
}





weights.r <- 1
rho.r <- 1

t <- as.matrix(traits.r)
# s <- as.matrix(traits.c)

# Some more abbreviations:
n.r <- dim(t)[1]
# n.c <- dim(s)[1]
# ntraits.r <- dim(t)[2]
# ntraits.c <- dim(s)[2]

# Make sure weights and abundances are of the right form, and
# introduce abbreviations g and h in line with theory:
g <- array(as.vector(weights.r),n.r)
# h <- array(as.vector(weights.c),n.c)
rho.r <- array(as.vector(rho.r),n.r)
# rho.c <- array(as.vector(rho.c),n.c)

# Normal approximation for resource trait distribution:
C.u <- cov.wt(t, wt=as.vector(rho.r*g), method="ML")
# (Note: Since we are using weights, "ML" is the better choice!)
mu.u <- C.u$center
C.u <- C.u$cov


# Check if standardization works for traits form crazy distributions and for the mixture
# of distributions

spec = 10
trai = 6

rbinom(trai, 1, p=0.7)
rnorm(trai,10, 2)
rexp(trai, 0.2)
runif(trai)

traitmat = matrix(0,nrow=spec, ncol = trai)
traitmat[,1]
