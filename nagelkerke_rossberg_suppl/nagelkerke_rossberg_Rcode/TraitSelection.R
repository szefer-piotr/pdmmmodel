## Find the most predictive combinations of phenotypic traits of
## consumers and resources.  This code takes a while, but can be run
## in batch mode.  All relevant output is written into a log file.

source("nichespace.R")
source("functions.R")

target.dim <- 4 # dimensionality of niche space

log.file.name <- paste("log", as.double(target.dim),".log",sep="")
log.file <- file(log.file.name, open = "w")

# Log dimensionality:
write(target.dim,file=log.file)
flush(log.file)

# Now log file is set up
#########################################################

# Read traits:
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

# Now we have all data together
###########################################################

# Parameters for cross.validate.a(..)
n.resources <- dim(resource.traits)[1]
n.consumers <- dim(consumer.traits)[1]
n.omitted.resources <- 1
n.omitted.consumers <- 1
n.reps <- n.resources*n.consumers
weights.r <- rep(1,n.resources)

# Exhaustive search of trait combinations for best cross.validate.a
best.cv <- 0
max.num.traits <- 2 # Max number of phenotypic traits per interaction partner 
for(rt in 1:max.num.traits){
  for(ct in 1:max.num.traits){
    for(c.traits in as.data.frame(combn(1:dim(consumer.traits)[2],ct))){
      print(c.traits) # phenotypic consumer traits to consider
      write(c.traits,file=log.file)
      for(r.traits in as.data.frame(combn(1:dim(resource.traits)[2],rt))){
        # combination of phenotypic traits to consider: 
        print(c(colnames(resource.traits)[r.traits],
                colnames(consumer.traits)[c.traits] ))
        # values of selected phenotypic traits:
        selected.resource.traits <-
          col.select(resource.traits,r.traits)
        selected.consumer.traits <-
          log10(col.select(consumer.traits,c.traits))
        # do cross-valication:
        this.cv <- cross.validate.a(dims.to.keep=target.dim)
        # log result:
        cat(c("XM ",this.cv,
              colnames(resource.traits)[r.traits],
              colnames(consumer.traits)[c.traits]),
            "\n",file=log.file)
        # if this is new optimum...
        if(this.cv > best.cv){
          # ... save and log new optimum:
          best.cv <- this.cv
          best.traits <-  list(r.traits=r.traits,c.traits=c.traits)
          cat("XB ",file=log.file)
          write(best.cv,file=log.file)
          cat(colnames(resource.traits)[r.traits],"\n",
              file=log.file)
          cat(colnames(consumer.traits)[c.traits],"\n",
              file=log.file)
          cat(r.traits,"\n",file=log.file)
          cat(c.traits,"\n",file=log.file)
        }
        print(this.cv)
      }
      print(best.traits)
      print(best.cv)
      print("")
    }
  }
}
quit(save = "no", runLast = FALSE)

# Note: the last line of the log file starting with "XB" contains the
# optimal trait combination
