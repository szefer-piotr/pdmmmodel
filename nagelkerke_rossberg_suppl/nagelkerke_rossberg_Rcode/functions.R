# $Id: functions.R,v 1.24.1.12 2010/10/03 13:47:51 axel Exp axel $	

## This selects colums of a matrix, and keeps the matrix structure even
## if only a single column is selected:
col.select <- function(m,cols){
  m1 <- as.matrix(as.data.frame(m)[,cols]);
  if(is.numeric(cols)){
    colnames(m1)<-colnames(m)[cols];
  }else{
    colnames(m1)<-cols;
  }
  rownames(m1)<-rownames(m);
  return(m1);
}

## This selects rows of a matrix, and keeps the matrix structure even
## if only a single row is selected:
row.select <- function(m,rows){
  m1 <- as.matrix(as.data.frame(m)[rows,]);
  if(is.numeric(rows)){
    rownames(m1)<-rownames(m)[rows];
  }else{
    rownames(m1)<-rows;
  }
  colnames(m1)<-colnames(m);
  return(m1);
}

## This generates a diagonal matrix, even for 1x1 matrices
my.diag <- function(v){
  if(length(v)==1){
    return(as.matrix(v))
  }else{
    return(diag(v))
  }
}

## Compute normalization of phenotypic traits to unit covariance.  The
## output are the normalized traits.  Important information for the
## back-transformation is saved in attributes of output data ("eigen",
## "traitnames", "mean.traits").
cov.normalize <- function(data, tol=1e-7, keep=ncol(data)){
  mean.traits <- colMeans(data)
  data <- t(t(data) - mean.traits)
#  print(data)
  c <- cov(data)
  e <- eigen(c)
  d <- data %*% e$vectors
  max.i <- ncol(data)
  for(i in 1:max.i){
    d[,i] <- d[,i]/sqrt(e$values[i])
    if(abs(e$values[i])<tol && i<=keep){
      warning("ill-conditioned covariance matrix")
    }
  }
  d <- d[,1:max.i]
  e$values <- e$values[1:max.i]
  e$vectors <- e$vectors[,1:max.i]
  attr(d,"eigen") <- e
  attr(d,"traitnames") <- colnames(data)
  attr(d,"mean.traits") <- mean.traits
  dim(d) <- dim(data)
  return(d)
}

## Randomly permute rows of a matrix
permute.species <- function(data){
  nspecies <- dim(data)[1]
  return(data[sample(1:nspecies,nspecies),])
}

## Compute weights at which different phenotypic traits contribute to
## trophic traits (diagnostic).
trophic.weigth.matrix <- function(trait.projector,e){
  M1 <- trait.projector %*% t(trait.projector)
  M2 <- e$vectors %*% my.diag(e$values) %*% t(e$vectors)
  w <- M1*M2
  return(w/sum(w))
}

## Extract trophic traits from fitted quadratic form and corresponding
## phenotypic traits.  Parameter sel is either an integer, indicating
## how many dimensions to keep, or (less useful) a list telling which
## dimensions to keep.
trophic.traits <- function(r.traits,c.traits,qform,sel=0,tol=1e-7){
  e <- eigen(qform$C) # eigenvectors and -values of C
  ##print(e)
  ## Some house keeping:
  r.dim <- dim(r.traits)[2]
  c.dim <- dim(c.traits)[2]
  n.r <- dim(r.traits)[1]
  n.c <- dim(c.traits)[1]
  v <- e$values # eigenvalues of C
  ## Get list of selected dimensions
  if(as.integer(sel)==sel){
    if(sel>0){
      sva <- sort.int(-abs(v),index.return=TRUE)
#      sva=sort.int(v,index.return=TRUE)
      sel <- sva$ix[1:sel]
      cat("selected indices ",sel,"\n")
    }else{
      sel=1:length(e$values)
    }
  }
  ## ... and its complement
  unsel <- setdiff(seq(length(e$values)),sel)
  ## Some housekeeping:
  dimension.names <- paste("dim", as.double(1:length(sel)),sep="")
  F.val <- array(0,dim=c(n.c,length(sel)))
  colnames(F.val) <- dimension.names
  V.val <- array(0,dim=c(n.r,length(sel)))
  colnames(V.val) <- dimension.names
  F.star <- vector(mode="numeric",length=n.c)*0
  V.star <- vector(mode="numeric",length=n.r)*0
  k <- 1
  a <- qform$a
#
  lambda.val <- e$value[sel]
    
  ### Compute "projectors", i.e. the linear transformations that lead
  ### from phenotypic traits to trophic traits:

  d <- attributes(r.traits)$eigen$values
  d <- 1/sqrt(abs(d))

  ## Compute projector onto vulnerability traits:
  V <- attr(r.traits,"eigen")$vectors %*%
    my.diag(d) %*%
      e$vectors[1:r.dim,sel] %*% my.diag(sqrt(abs(lambda.val)))
  rownames(V) <- attr(r.traits,"traitnames")
  V.offset <- -attr(r.traits,"mean.traits") %*% V
  colnames(V.offset) <- dimension.names

  
  ## Canonicalize signs of trait vectors
  if(V[1,1]<0){
    e$vectors <- -e$vectors
    V <- -V
    V.offset <- - V.offset
  }

  ### Now the same for foraging traits
  
  d <- attributes(c.traits)$eigen$values
  d <- 1/sqrt(abs(d))
  
  ## Compute projector onto foraging traits:
  F <- -attr(c.traits,"eigen")$vectors %*%
    my.diag(d) %*%
      e$vectors[r.dim+(1:c.dim),sel] %*% my.diag(sqrt(abs(lambda.val)))
  rownames(F) <- attr(c.traits,"traitnames")
  F.offset <- -attr(c.traits,"mean.traits") %*% F
  colnames(F.offset) <- dimension.names
  
  for(i in sel){
    ## Compute trophic traits and report results:

    print(paste("Dimension ", as.double(k),":",sep=""),quote=FALSE)
    print(paste("lambda ", as.double(e$value[i]),":",sep=""),quote=FALSE)
    print("V:",quote=FALSE)
    print(t(V[,k]))

    print("F:",quote=FALSE)
    print(t(F[,k]))
    
    V.val[,k] <- r.traits %*% e$vectors[1:r.dim,i]
    F.val[,k] <- -c.traits %*% e$vectors[r.dim+(1:c.dim),i]
    V.val[,k] <- V.val[,k] +
      niche.dot(qform$b[],e$vectors[,i])/lambda.val[k]
    V.offset[k] <- V.offset[k] +
      niche.dot(qform$b[],e$vectors[,i])*sign(lambda.val[k])/
        sqrt(abs(lambda.val[k]))
    V.val[,k] <- V.val[,k] * sqrt(abs(lambda.val[k]))
    F.val[,k] <- F.val[,k] * sqrt(abs(lambda.val[k]))

    a <- a - niche.dot(qform$b,e$vectors[,i])^2/(2*lambda.val[k])
    
    k <- k+1
  }

  for(i in unsel){
    weight <- niche.dot(qform$b,e$vectors[,i])
    V.star <- V.star + r.traits %*% (weight*e$vectors[(1:r.dim),i])
    F.star <- F.star + c.traits %*% (weight*e$vectors[r.dim+(1:c.dim),i])
  }

  colnames(V)<-colnames(V.val)<-dimension.names
  colnames(F)<-colnames(F.val)<-dimension.names

  return(list(V.projector=V,V.offset=V.offset,
              F.projector=F,F.offset=F.offset,
              V=V.val,F=F.val,
              V.star=V.star,F.star=F.star,
              lambda=lambda.val, sigma=-sign(lambda.val), a=a ))
}

## Helper function to delete effect of selected niche-space dimensions
## from quadratic form.
delete.dims <- function(sel,qf){
  e <- eigen(qf$C)
  v <- e$values
  if(length(sel)==1 && is.integer(sel)){
    if(sel>0){
      sva=sort.int(-abs(v),index.return=TRUE)
      v[sva$ix[-(1:sel)]] <- 0
    }else{
      v[] <- 0
    }
  }else if(is.vector(sel)){
    v[-sel] <- 0
  }
  qf$C=e$vectors %*% my.diag(v) %*% t(e$vectors)
  return(qf)
}


## This does the hard work of cross-validation.  This is done taking
## some short cuts, i.e. without actually computing the trophic
## traits.  Some global variables have an effect here, among these:
## selected.resource.traits, selected.consumer.traits, n.reps,
## n.omitted.consumers, n.omitted.resources,
cross.validate.a <- function(dimdeletion=NULL,dims.to.keep=NULL,
                             r.comp,c.comp,tol=1e-7){

  ## We start with lots of book keeping:
  if(missing(r.comp)){
    r.comp <- ncol(selected.resource.traits)
  }
  if(missing(c.comp)){
    c.comp <- ncol(selected.consumer.traits)
  }
  n.resources <- dim(selected.resource.traits)[1]
  n.consumers <- dim(selected.consumer.traits)[1]
  a.corr <- vector("numeric",
                   length=n.reps
                   )
  all.aa <- c();
  all.aa.true <- c();
  all.wt <- c();
  k <- 1
  for(rep in 1:n.reps){
    if(n.reps==n.consumers*n.resources &&
       n.omitted.consumers==1 &&
       n.omitted.resources==1
       ){
      omitted.resources <- ((rep-1)%/%n.consumers)+1
      omitted.consumers <- ((rep-1)%%n.consumers) +1
    }else
    if(n.reps==n.resources &&
       n.omitted.consumers==0 &&
       n.omitted.resources==1
       ){
      omitted.consumers <- NULL
      omitted.resources <- rep:rep
    }else{
      omitted.resources <- sample(1:n.resources,n.omitted.resources)
      omitted.consumers <- sample(1:n.consumers,n.omitted.consumers)
    }
    kept.resources <- setdiff(1:n.resources,omitted.resources)
    kept.consumers <- setdiff(1:n.consumers,omitted.consumers)
                                        #
    input.resource.traits <-
      cov.normalize(row.select(selected.resource.traits,kept.resources),keep=r.comp)
    input.consumer.traits <-
      cov.normalize(row.select(selected.consumer.traits,kept.consumers),keep=c.comp)

    ## Here the actual work start.  The quadratic form is computed
    ## based on the select data.
    QF <- niche.quadratic.form(input.resource.traits,
                               input.consumer.traits,
                               interactions[kept.resources,kept.consumers]
                               ,weights.r=weights.r[kept.resources]
                               ,rho.r=res.abundance[kept.resources]
                               )
    ## Remove information for dimensions that are not included:
    if(length(dimdeletion)>0){
      QF <- delete.dims(as.vector(dimdeletion),QF)
    }
    if(length(dims.to.keep)>0){
      QF <- delete.dims(as.integer(dims.to.keep),QF)
    }
                                        #

    d <- attributes(input.resource.traits)$eigen$values
    d <- 1/sqrt(abs(d))
    r.normalizer <-
      attributes(input.resource.traits)$eigen$vectors %*% my.diag(d)
    
    d <- attributes(input.consumer.traits)$eigen$values
    d <- 1/sqrt(abs(d))
    c.normalizer <-
      attributes(input.consumer.traits)$eigen$vectors %*% my.diag(d)

    # Now use the omitted consumers and resources (usually just one
    # each) to see how good the model predicts their interaction
    # strengths:
    for(r in omitted.resources){
      for(c in omitted.consumers){
        aa <- niche.a((selected.resource.traits[r,] -
                       attributes(input.resource.traits)$mean.traits ) %*%
                      r.normalizer,
                      (selected.consumer.traits[c,] -
                       attributes(input.consumer.traits)$mean.traits ) %*%
                      c.normalizer,
                      QF)
        aa.true <- interactions[r,c]/res.abundance[,][r]

        if(is.finite(aa) && is.finite(aa.true)){
          k <- k+1;
          all.aa <- c(all.aa,aa)
          all.aa.true <- c(all.aa.true,aa.true)
          all.wt <- c(all.wt,res.abundance[,][r]^2)
        }
      }
    }
  }
  return(cov.wt(cbind(all.aa,all.aa.true),all.wt,cor = TRUE)$cor[1,2])
}
