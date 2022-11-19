
########################## define functions for simulation study ####################################

# outliers sample function
outl.sample <- function(p, N, corr.eps, gamma, F.mat){
  
  # constants
  quant = 0.95
  
  # retype
  F.mat <- as(F.mat,"matrix")
  
  # V.mat
  V.mats      = V.mat.fun(F.mat)
  V.mat       = V.mats$V.mat
  V2.mat.proj = V.mat %*% t(V.mat)
  constrs     = V.mats$V.mat.o
  dim.ker     = p-dim(V.mat)[2]
  
  
  # covariance matrix 
  inds    = combn(1:p,2)
  cov.mat = matrix(0, p, p)
  diag(cov.mat) = 1
  cov.mat[lower.tri(cov.mat)] = (-0.9)^abs(inds[1,]-inds[2,])
  cov.mat = 0.5 * (cov.mat + t(cov.mat))
  
  # transformed covariance matrix
  cov.mat.trans     = V.mat %*% t(V.mat) %*% (cov.mat) %*% V.mat %*% t(V.mat)
  cov.mat.trans.inv = V.mat %*% t(V.mat) %*% ginv(cov.mat) %*% V.mat %*% t(V.mat)
  #cov.mat.trans = ginv(cov.mat.trans)
  
  # generate samples
  x.mat = mvrnorm(N, rep(0,p), cov.mat.trans)
  
  ### corrupt data 
  # generate indices - first column p.idx, second N.idx
  ind.corr = Reduce("rbind",lapply(1:p, function(p.idx){
    cbind(sample(1:N, ceiling(N*corr.eps), replace = F), rep(p.idx,ceiling(N*corr.eps)))}))
  colnames(ind.corr) = c("n.idx", "p.idx")
  ind.corr = ind.corr[order(ind.corr[,1], decreasing = F),]
  
  # replace some with corrupted samples
  x.corr.mat = lapply(1:N, function(n.idx){ 
    
    ind.rep.loc  = ind.corr[which(ind.corr[,"n.idx"] == n.idx),"p.idx"]   # the indices for given n.idx to corrupt 
    if(length(ind.rep.loc) != 0){
      
      # indices to contaminate
      k = length(ind.rep.loc)
      
      # generate u
      u.vec              = rep(0,p)
      u.vec[ind.rep.loc] = runif(length(ind.rep.loc),-1,1)
      
      # compute ---> 1/mu.resc^2 * ac + 1/mu.resc * bc + cc = 0
      ac = t(u.vec) %*% (cov.mat.trans.inv %*% u.vec)
      bc = -2 * t(u.vec) %*% (cov.mat.trans.inv %*% x.mat[n.idx,])
      cc = -qchisq(quant,p - dim.ker) + t(x.mat[n.idx,]) %*% cov.mat.trans.inv %*% x.mat[n.idx,]
      
      if((bc^2-4*ac*cc) < 0){
            return(x.mat[n.idx,])
      }else{
            # rescale so that  MD = qchisq(quant, p-1)
            mu.resc = as.numeric(2*ac/(-bc - sqrt(bc^2-4*ac*cc)))
            u.vec = u.vec / mu.resc
            x.loc = x.mat[n.idx,] - gamma * u.vec
            # t(x.mat[n.idx,]-u.vec) %*% cov.mat.trans.inv %*% (x.mat[n.idx,]-u.vec) # qchisq(quant, p-1)
      }
      
      return(x.loc)
    }else{
      return(x.mat[n.idx,])
    }
  })
  x.corr.mat = Reduce("rbind", x.corr.mat)
  
  # check which ones are  
  #apply(x.corr.mat, 1, function(u){ t(u)%*% ginv(cov.mat.trans) %*% u }) > qchisq(0.95,p-1)
  
  
   # add noise in orthogonal direction 
  ortho.mult      = rnorm(N, 0, 10)
  x.corr.mat.orig = x.corr.mat 
  x.mat.orig      = x.mat
  x.diff.mat = x.corr.mat - x.mat
  x.corr.mat = x.corr.mat %*% V2.mat.proj + ortho.mult %*% t(constrs) 
  x.mat      = x.mat %*% V2.mat.proj #+ ortho.mult %*% t(constrs) 

  # transform index of outliers to vector  
  ind.corr =  sort((ind.corr[,2]-1)*N + ind.corr[,1])
  
  return(list(x.mat = x.mat, 
              x.corr.mat.orig = x.corr.mat.orig,
              x.mat.orig      = x.mat.orig,
              x.mat           = x.mat,
              x.corr.mat = x.corr.mat, 
              x.diff.mat = x.diff.mat,
              ind.corr = ind.corr, 
              cov.mat.trans = cov.mat.trans))
}


PRFscores <- function(est.inds, real.inds){
  
  inter.len = length(intersect(est.inds,real.inds)) 
  real.len  =  length(real.inds)
  est.len   = length(est.inds)
  
  # precision
  if(real.len!=0){ precision = inter.len / real.len }else{ precision = 0}
  
  #recall
  if(est.len!=0){ recall = inter.len / est.len }else{ recall = 0}

  # F-score
  if(precision != 0 | recall != 0){
    Fscore = 2 * (precision * recall) / (precision + recall)
  }else{
    Fscore = 0
  }
  
  return(unlist(list(recall = recall, precision = precision, Fscore = Fscore)))
}




