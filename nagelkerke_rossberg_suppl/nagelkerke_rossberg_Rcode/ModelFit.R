## Computations to construct trophic niche space for given selection
## of phenotypic traits.

source("nichespace.R")
source("functions.R")

## Read traits:
resource.traits <-
  as.matrix(read.table("Tana_Resource_Traits.dat",header=TRUE))
consumer.traits <-
  as.matrix(read.table("Tana_Consumer_Traits.dat",header=TRUE))

# Read interaction strength

# Note: these are interaction per consumer individual, so we have to
# assume consumer densities = 1

diet.fractions <- read.table("Tana_Interaction_Strength.dat")/100

interactions <- t(diet.fractions)

# Guess relative resources abundance from stomach content
res.abundance <- t(t(colMeans(as.matrix(diet.fractions))))

## Now we have all data together
###########################################################

#### Parameters for cross.validate.a(..)
## number of resources
n.resources <- dim(resource.traits)[1] 
## number of consumers
n.consumers <- dim(consumer.traits)[1]
## how many to omit in cross-validation
n.omitted.resources <- 1
n.omitted.consumers <- 1
## how many repetitions (matters when using randomization)
n.reps <- n.resources*n.consumers
## option to assign different statistical weights to resources
weights.r <- rep(1,n.resources)


## Specify and double-check a set of dims + expanatory variables
## Specify by name:
source("functions.R")
dims.to.keep <- 4
r.list <- "R_pelagic R_macrodim"
c.list <- "Barbel Fork_length"
selected.resource.traits <-
  col.select(resource.traits,unlist(strsplit(r.list," ")))
selected.consumer.traits <-
  log10(col.select(consumer.traits,unlist(strsplit(c.list," "))))
## Verify cross valiation:
cross.validate.a(dims.to.keep=dims.to.keep)

# Fit the model
source("nichespace.R")
input.resource.traits <-
  cov.normalize(selected.resource.traits)
input.consumer.traits <-
  cov.normalize(selected.consumer.traits)
QF0 <- niche.quadratic.form(input.resource.traits,
                           input.consumer.traits,
                           interactions
                           ,weights.r=weights.r
                           ,rho.r=res.abundance[,]
                           )
print(eigen(QF0$C)$values)

# Truncate niche-space dimensions:
QF <- delete.dims(as.integer(dims.to.keep),QF0)

# Compute predicted link strengths:
link <- vector("character",
               length=dim(input.consumer.traits)[1]*
               dim(input.resource.traits)[1]
               ) # descriptions of links
aa <- vector("numeric",
             length=dim(input.consumer.traits)[1]*
             dim(input.resource.traits)[1]
             ) # list of predicted link strengths
aa.true <- aa; # list of observed link strengths
AA <- t(diet.fractions*0) # matrix of link strengths
predicted.diet <- array(0,dim=c(n.consumers,n.resources))
k <- 1;
for(c in 1:dim(input.consumer.traits)[1]){
  for(r in 1:dim(input.resource.traits)[1]){
    a <- niche.a(input.resource.traits[r,],
                 input.consumer.traits[c,],
                 QF)
    AA[r,c] <- a
    aa[k] <- a
    aa.true[k] <- interactions[r,c]/res.abundance[,][r]
    predicted.diet[c,r] <-  aa[k]*res.abundance[,][r]
    link[k]=paste(as.double(r), "->", as.double(c),sep="")
    k <- k+1;
  }
  predicted.diet[c,] <- predicted.diet[c,]/sum(predicted.diet[c,])
}

# Compute weighted correlation:
cov.wt(cbind(as.vector(as.matrix(aa)),as.vector(aa.true)),wt=rep(res.abundance,n.consumers)^2,cor=TRUE)$cor[1,2]

# Compute plot of predicted vs observed diets, and R^2:
plot(as.vector(t(t(diet.fractions[,]))),as.vector(predicted.diet[,]),xlim=c(0,1),ylim=c(0,1))
R2 <- cor(as.vector(t(t(diet.fractions[,]))),as.vector(predicted.diet[,]))^2
R2

# Same plot on log axes:
plot(as.vector(t(t(diet.fractions[,]))),as.vector(predicted.diet[,]),log="xy",xlim=c(0.001,1),ylim=c(0.001,1))

# Re-construct trophic trait space:
t.traits <- trophic.traits(input.resource.traits,input.consumer.traits,QF0,sel=dims.to.keep)

# The list t.traits now contains all info to work with and interpret
# trophic traits.
#############################################################

## Verify trophic distance matrix
distances <- matrix(nrow=nrow(t.traits$V),ncol=nrow(t.traits$F))
V <- t.traits$V  # not acutally permuted!!
F <- t.traits$F  # not acutally permuted!!
val.list <- matrix(nrow=n.resources*n.consumers,ncol=6)
for(r in seq(nrow(V))){
  for(c in seq(nrow(F))){
    ddd <- V[r,]-F[c,]
    distances[r,c]=-sum(ddd^2*sign(t.traits$lambda))-
      2*(t.traits$V.star[r]+t.traits$F.star[c])
  }
}
a.predicted <- exp(t.traits$a-distances/2)
plot(log10(as.vector(as.matrix(a.predicted))),log10(as.vector(aa)))
