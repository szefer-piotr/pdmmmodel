# Statistical propoerties of link strength matrix
source("nagelkerke_rossberg_suppl/nagelkerke_rossberg_Rcode/nichespace.R")
source("nagelkerke_rossberg_suppl/nagelkerke_rossberg_Rcode/functions.R")
n.r = 25
n.c = 36
ntraits.r=10
ntraits.c=10

# Random normal trait data
mock <- niche.mock.data(n.r = 25,
                        n.c = 36,
                        ntraits.r,
                        ntraits.c,
                        b=0,
                        C=NA)

v <- mock$traits.r
# Normalize covariance
v <- cov.normalize(v)
f <- mock$traits.c

Qv <- cov(v)
Qf <- cov(f)

# Covariance matrix for differences v - f for a random pair of species
Q <- Qv + Qf
S <- diag(x = sample(c(-1,1), size = 10, replace = T), 10,10) # signature

