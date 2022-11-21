library(limSolve)
library(genlasso)
library(ICSNP)
library(MASS)
library(lars)
library(cellWise)
library(gglasso)
library(robustbase)





# function to calculate the graph difference matrix for given 
# p ... original dimension
# ed.mat...edges (2*N matrix) 
# ws.vec ...weights (1*N vector)
diff.mat.fun <- function(p, ed.mat, ws.vec){
  
  # init 
  edge.nbr = dim(ed.mat)[2]
  diff.mat.loc = Matrix(0, nrow =edge.nbr, ncol = p)
  
  # set weights in diff.matrix
  diff.mat.loc[cbind(1:edge.nbr,ed.mat[1,])] = ws.vec
  diff.mat.loc[cbind(1:edge.nbr,ed.mat[2,])] = -ws.vec

  return(diff.mat.loc)
}


# svd of F.mat
F.mat.svd <- function(F.mat){
  
  # compute the eigenvecs
  tol = sqrt(.Machine$double.eps)
  D.svd    = svd(x = F.mat, nv = min(dim(F.mat)), nu = 0)
  positive = min(sum(D.svd$d >  max(tol * D.svd$d, 0)), min(dim(F.mat)))

  return(list(D.mat = D.svd$d[1:positive], V.mat = as.matrix(D.svd$v[,1:positive])))
}



# function to get the coordinate map (clr but for the F.mat matrix)
coord.map <- function(x.mat, V.mat){
  
  # transform
  x.trans = x.mat %*% (V.mat %*% t(V.mat))
  
  # names
  colnames(x.trans) = colnames(x.mat)
  rownames(x.trans) = rownames(x.mat)
  
  return(x.trans)
}


# rescaling function
rescaling <- function(x.mat, F.svd){
  
  # init
  const.med = 1.4826
  tol = sqrt(.Machine$double.eps)
  N  = dim(x.mat)[1]
  px = dim(x.mat)[2]
  pD = dim(F.mat)[1]
  
  # compute one to one transform
  x.trans.svd = x.mat %*% (F.svd$V.mat %*% diag(F.svd$D.mat))
   
  # compute robust center 
  x.centers = spatial.median(x.trans.svd) # in transformed coordinates
  x.centers = F.svd$V.mat %*% (diag(1/F.svd$D.mat) %*% x.centers)
  
  # center data
  x.mat.new = scale(x.mat, center = x.centers, scale = F) # sapply(1:px, function(p.idx){ x.mat[,p.idx] - x.centers[p.idx]})
  colnames(x.mat.new) = colnames(x.mat)
  rownames(x.mat.new) = rownames(x.mat)
  
  # now compute scales
  Dx.scales = sapply(1:pD, function(pd.idx){ const.med * median(abs(F.mat[pd.idx,] %*% t(x.mat.new)))  }) 
  F.mat.new = diag(1 / Dx.scales) %*% F.mat
  
  
  return(list(x.mat.new = x.mat.new, F.mat.new = F.mat.new, x.centers = x.centers, Dx.scales = Dx.scales))
}


# calculate the precision matrix from from x and V.mat 
P.mat.fun <- function(x.mat, col.centers, V.mat){
  
  x.centered  = scale(x.mat, center = col.centers, scale = F)
  x.trans     = t(t(V.mat) %*%  t(x.centered))
  x.trans.inv = red.ginv(1 / (dim(x.trans)[1]-1) * t(x.trans) %*% x.trans )
  P.mat.sq    = V.mat %*% x.trans.inv$InvSqrt %*% t(V.mat)
  P.mat       = V.mat %*% x.trans.inv$Inv %*% t(V.mat)
  S.mat       = V.mat %*% cov(x.trans) %*% t(V.mat)
  

  return(list(P.mat = P.mat, P.mat.sq = P.mat.sq, S.mat = S.mat))
}

# function to calculate generalized inverse of matrix 
red.ginv <- function(A){
  
   # invert
   tol = sqrt(.Machine$double.eps)
   x.trans.svd = svd(x = A, nv =  min(dim(A))) 
   positive     = 1:min(max(which(x.trans.svd$d >  max(tol * x.trans.svd$d, 0))))
   #positive    = x.trans.svd$d >  max(tol * x.trans.svd$d, 0)
   A.inv       = x.trans.svd$v[,positive, drop = FALSE] %*% diag(1 / x.trans.svd$d[positive, drop = FALSE]) %*% t(x.trans.svd$v[,positive, drop = FALSE]) 
   A.inv.sqrt  = x.trans.svd$v[,positive, drop = FALSE] %*% diag(1 / x.trans.svd$d[positive, drop = FALSE])^(1/2) %*% t(x.trans.svd$v[,positive, drop = FALSE]) 

  return(list(Inv = A.inv, InvSqrt = A.inv.sqrt))
}



#least squares
gen.ols <- function(x, y, idx, eps){
 
  # get reduced constraints
  loc.coeffs = lm.fit(x = rbind(matrix(x[,idx], ncol = length(idx)), sqrt(eps) * diag(rep(1,length(idx)))), 
                          y = c(y, rep(0,length(idx))))$coefficients

  
  return(loc.coeffs) 
}

# lasso with groups
lasso.fun <- function(predictors, response, weights.loc, grp.df){
  
  # set some coeffs to zero that are too small
  eps.zero  = 1e-8
  eps.ols   = 1e-2
  
  # init 
  p      = dim(predictors)[2]
  groups = unique(grp.df$groups)
  n.grps = length(groups)
  
  ### solve with gglasso 
  # copy dataframe columns for same group
  ord = order(grp.df$groups)
  groupord = grp.df$groups[ord]
  varord   =  grp.df$vars[ord]
  predictors.ext = predictors[, varord]

  # run glasso on copied data 
  weights.loc.inv = 1 / weights.loc[grp.df$groups]
  mod = gglasso(x = predictors.ext %*% diag(weights.loc.inv), 
                y     = response, 
                group = groupord, 
                loss  = "ls",
                intercept = F,
                nlambda   = 400)  
  mod$beta = diag(weights.loc.inv) %*% mod$beta
  mod$var  =  apply(mod$beta, 2, FUN = function(x) {
        varord[which(x != 0)] })

  
  # reduce to get coefficient path - groups that enter
  nnzero.path = sapply(1:length(mod$var), function(l.idx){ 
    if(l.idx == 1){ return(NULL) }
    val = union(mod$var[[l.idx]],mod$var[[l.idx-1]])
    val = setdiff(val, mod$var[[l.idx-1]])
    if(length(val)==0){return(NA)}else{ return(val) }
  }, simplify = F)
  nnzero.path = nnzero.path[!is.na(nnzero.path)]
  nnzero.path = Reduce(f = function(u,v){ union(u,v) }, nnzero.path, accumulate = TRUE)
  # nnzero.path = c(nnzero.path, list(1:p))  # adding the full solution as well
  ordering    = Reduce(f = function(u,v){ unique(c(u,v)) }, nnzero.path)

  addedvars   = sapply(2:length(nnzero.path), function(l.idx){
    setdiff(nnzero.path[[l.idx]], nnzero.path[[l.idx-1]]) 
  }, simplify  = F)
  
    
  # assign final beta - get OLS solutions 
  beta = lapply(nnzero.path, function(idx){ 
    
    if(length(idx)==0){ 
      return(rep(0,p)) 
    }else{
      
      # solve OLS
      coeffs     = rep(0,p)
      loc.coeffs = gen.ols(predictors, response, idx, eps.ols)  
      coeffs[idx]= loc.coeffs
      
      return(coeffs)
    }
    
    })
  beta = as.matrix(Reduce("rbind", beta))

  #  get RSS
  RSS = apply(beta, 1, function(u){ sum((response - predictors %*% u)^2) })
  
  return(list(RSS = RSS, 
              beta = beta,
              ordering = ordering, 
              nnzero.path = nnzero.path,
              addedvars = addedvars))
}


# groups is a two col dataframe with first col the variable index and second
# col the group it belongs to 
cellHandlerExt <- function(X, F.mat, grp.df, quant = 0.99){
  
  # TODO: implement check of extension
  # badinds = which(abs(sdres) > sqrt(qchisq(quant, 1)))   
  # Todo: groups must be consecutive! without any jumps from 1 to 5 and no groups 2,3,4
  
  
  # reorder groupings in ascending order
  ord      = order(grp.df$groups)
  groupord = grp.df$groups[ord]
  varord   = grp.df$vars[ord]
  grp.df$vars   = varord
  grp.df$groups = groupord
  
  # dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  
  # groups infos 
  groups   = unique(grp.df$groups)
  n.groups = length(groups)
  
  
  # compute svd 
  F.svd = F.mat.svd(F.mat)
  V.mat = F.svd$V.mat
  F.rk  = dim(V.mat)[2]
  
  # compute rescaling etc
  res.obj   = rescaling(X, F.svd)
  x.centers = res.obj$x.centers
  X.new     = res.obj$x.mat.new
  F.new     = res.obj$F.mat.new
  
  # 
  V2.mat = (F.svd$V.mat %*% t(F.svd$V.mat))
  X.proj = X.new %*% (diag(rep(1,p))-V2.mat)
  X.new  = X.new %*% V2.mat

  # compute initial mean/cov estimates
  init.sol = cellWise:::DDCWcov(X.new)$Z
  mu       = colMeans(init.sol)
  P.out    = P.mat.fun(init.sol, mu, V.mat)
  Sigma    = P.out$S.mat
  Prec     = P.out$P.mat     
  Prec.sq  = P.out$P.mat.sq
  predictors = Prec.sq

  # compute the precision matrices for each grouping 
  prec.groups = lapply(seq_len(n.groups), function(grp.idx){
     var.idx = grp.df$vars[grp.df$groups == groups[grp.idx]]
     ginv(Sigma[var.idx,var.idx])
  })
  
  # compute the weight estimates - for each group
  weights = sapply(seq_len(n.groups), function(grp.idx){
    var.idx = grp.df$vars[grp.df$groups == groups[grp.idx]]
    xs      = scale(X.new[,var.idx], center = mu[var.idx], scale = F)
    
    # mahalanobis distance
    mhd = sqrt(diag(xs %*% prec.groups[[grp.idx]] %*%  t(xs)))
    1.5 / mhd # pmin(1,1.5 / mhd)
  }, simplify = "matrix")
  

   # init next 
   Ximp     = X.new
   Zres     = matrix(0, n, p)
   Zres.den = matrix(0, n, p)
   Zres.num = Zres.den
   # run detection one time 
   for (i in seq_len(n)){
        x            = X.new[i, ]
        response     = predictors %*% (x - mu)
        weights.loc  = weights[i,]
        
        # use groupwise lasso
        lassOut = lasso.fun(predictors, response, weights.loc, grp.df)
        cat("lass: ", lassOut$ordering, "\n")
        cat("lass: ", round(lassOut$RSS,2),"\n")
        
        #zz = lars.ext.unconstr(predictors, response, constrs, V2.mat.proj, weights.loc)
        #cat("lars:", zz$ordering, "\n")
        #cat("lars:", zz$RSS, "\n")
        
        # path index of possible outlying cells 
        # deltas   =  abs((diff(lassOut$RSS)))
        # badCells = which(deltas > qchisq(quant, 1))
        
        badCells = which(lassOut$RSS > qchisq(quant, F.rk))

        # check bad cells and impute 
        if(length(badCells) > 0){
            
            # corrupted cells indices
            ind     = min(max(badCells)+1,length(lassOut$nnzero.path))
            badinds = lassOut$nnzero.path[[ind]] # take at most the allowed one 
            
            
            
            # if(length(badinds) == p){
            #     cond.mhds = (x - mu) / sqrt(diag(Sigma))
            # }else{
            #    
            #   
            #    # compute conditional covariance 
            #    con.cov = Sigma[badinds, badinds] - Sigma[badinds,-badinds]%*% ginv(Sigma[-badinds,-badinds])%*%Sigma[-badinds,badinds]
            #    
            #    if(length(badinds) == 1 && abs(con.cov) < 1e-4){ 
            #      cond.mhds = rep(0, length(badinds.added)) 
            #     }
            #    
            #    # check for each  subgroups really outliers
            #    badinds.added = lassOut$addedvars[1:max(badCells)]
            # 
            #    # compute conditional mahalanobis distances
            #    cond.mhds = lapply(badinds.added, function(o.idx){ 
            #    
            #      ind      = which(badinds == o.idx)
            #      beta.loc = lassOut$beta[max(badCells)+1,o.idx]
            #      sqrt(t(beta.loc) %*% ginv(con.cov[ind,ind]) %*% (beta.loc))
            #         
            #   })
            #   cond.mhds = Reduce("c", cond.mhds) 
            #    
            # }
            # 
            # # recheck groups 
            # badinds = which(abs(sdres) > sqrt(qchisq(quant, 1)))   
            # 
          
            # impute data 
            if(length(badinds) == p){
                              
                Ximp[i, ]     = mu
                Zres.num[i, ] = (x - mu)
                Zres.den[i, ] = sqrt(diag(Sigma))
                Zres[i, ]     = Zres.num[i, ]/Zres.den[i, ]
            
            }else{
              
                sdres = Sigma[badinds,badinds] - Sigma[badinds,-badinds]%*% ginv(Sigma[-badinds,-badinds])%*%Sigma[-badinds,badinds]
                if(abs(sdres) < 1e-6 && length(badinds)==1){ next }
              
                # get best prediction
                best.guess = gen.ols(predictors, response, badinds, eps.ols)
                
                # impute
                ximp          = X.new[i,]
                ximp[badinds] = ximp[badinds] - best.guess 
                Ximp[i, ]     = ximp
                
                # calculate residuals and Zres.num/Zres.den
                residual             = X.new[i, ] - ximp
                Zres.num[i, badinds] = residual[badinds]
                Zres.den[i, badinds] = sqrt(diag(abs(sdres)))  
                Zres[i, badinds]     = Zres.num[i, badinds] / Zres.den[i,badinds]

            }
        }
    }
    
    # get outlying cells indices 
    indCells = which(abs(Zres) > sqrt(qchisq(quant, 1)))

    # names
    colnames(Zres) = colnames(X); rownames(Zres) = rownames(X)
    colnames(Ximp) = colnames(X); rownames(Ximp) = rownames(X)

    # add back projected part
    Ximp.orig = Ximp + X.proj
        
    return(list(Ximp = Ximp, Ximp.orig = Ximp.orig, indCells = indCells, 
        Zres = Zres, Zres.den = Zres.den))
}



  x.corr.mat = log(robCompositions::expendituresEU)
  x.corr.mat = as.matrix(x.corr.mat)
  p = dim(x.corr.mat)[2]
 
  F.mat = diff.mat.fun(p,combn(1:p,2),rep(1,p*(p-1)/2))

  grp.df = cbind(1:12,1:12)
  #grp.df = cbind(1:5,1:5)
  grp.df = data.frame(grp.df)
  names(grp.df) = c("vars","groups")
  

  bla = cellHandlerExt(x.corr.mat, F.mat, grp.df, 0.95)
  plot(abs(!!bla$Zres))
  plot(abs(bla$Zres))

  
  