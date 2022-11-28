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


# calculate betas from a list of variabels (nnzero.path)
beta.fun <- function(nnzero.path, x, mu, precision){
  
  if(!is.list(nnzero.path)){ nnzero.path = list(nnzero.path)}
  
  beta = lapply(nnzero.path, function(idx){ 
    
    if(length(idx)==0){ 
      return(rep(0,p)) 
    }else{
      
      # solve OLS
      badinds    = idx # the variables to replace
      coeffs     = rep(0,p)
      loc.coeffs = x[badinds]-mu[badinds]
      loc.coeffs = loc.coeffs + ginv(precision[badinds,badinds]) %*% precision[badinds,-badinds] %*% (x[-badinds]-mu[-badinds])
      coeffs[badinds]= loc.coeffs
      
      return(coeffs)
    }
    
    })
  if(length(nnzero.path) == 1){
    return(unlist(beta))
  }else{
    return(as.matrix(Reduce("rbind", beta))) 
  }
  

}


# function that takes a matrix of coefficients and calculates the RSS
RSS.fun <- function(beta, response, predictors){
  
  RSS = apply(beta, 1, function(u){ sum((response - predictors %*% u)^2) })

  return(RSS)
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
  
  # adding the full solution as well
  if(!(length(nnzero.path[[length(nnzero.path)]]) == p)){
    nnzero.path = c(nnzero.path, list(1:p))    
  }
  
  # reduce the solution path so that only unique ones are kept
  ordering = Reduce(f = function(u,v){ unique(c(u,v)) }, nnzero.path)

  # added variables in each step
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
singular.cellImputation <- function(X, V.mat, grp.df, mu, B.mat, alpha = 0.99){
  
  
  # TODO:  do subblock check
  
  # reorder groupings in ascending order
  # ord      = order(grp.df$groups)
  # groupord = grp.df$groups[ord]
  # varord   = grp.df$vars[ord]
  # grp.df$vars   = varord
  # grp.df$groups = groupord
  # 
  
  # dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  
  # groups infos 
  groups   = unique(grp.df$groups)
  n.groups = length(groups)
  
  
  # compute svd 
  # F.svd = F.mat.svd(F.mat)
  # V.mat = F.svd$V.mat
  F.rk  = dim(V.mat)[2]
  
  # calculate 
  V2.mat = V.mat %*% t(V.mat)
  X.proj = X %*% (diag(rep(1,p))-V2.mat)
  X.new  = X %*% V2.mat
  
  # compute rescaling etc
  # res.obj   = rescaling(X.new, F.svd)
  # x.centers = res.obj$x.centers
  # X.new     = res.obj$x.mat.new
  # F.new     = F.mat#res.obj$F.mat.new


  # compute initial mean/cov estimates
  Sigma      = V.mat %*% B.mat %*% t(V.mat)
  B.out      = red.ginv(B.mat)
  precision  = V.mat %*% B.out$Inv %*% t(V.mat)
  predictors = V.mat %*% B.out$InvSqrt %*% t(V.mat)
  

  # compute the precision matrices for each grouping 
  prec.groups = lapply(seq_len(n.groups), function(grp.idx){
     var.idx = grp.df$vars[grp.df$groups == groups[grp.idx]]
     precision[var.idx,var.idx]
  })
  
  # compute the weight estimates - for each group
  weights = sapply(seq_len(n.groups), function(grp.idx){
    var.idx = grp.df$vars[grp.df$groups == groups[grp.idx]]
    xs      = scale(X.new[,var.idx], center = mu[var.idx], scale = F)
    
    # mahalanobis distance
    mhd = sqrt(diag(xs %*% prec.groups[[grp.idx]] %*%  t(xs)))
    pmin(1,1.5 / mhd) #1.5 / mhd #
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
        lassOut$beta = beta.fun(lassOut$nnzero.path, x, mu, precision)
        lassOut$RSS  = RSS.fun(lassOut$beta, response, predictors)
        
        cat("lass: ", lassOut$ordering, "\n")
        cat("lass: ", round(lassOut$RSS,2),"\n")
        
        #lassOut = lars.ext.unconstr(predictors, response, V.mat %*% t(V.mat), weights.loc)
        #cat("lars:", zz$ordering, "\n")
        #cat("lars:", round(zz$RSS,2), "\n")
        
        # path index of possible outlying cells 
        #deltas   =  abs((diff(lassOut$RSS)))
        #badCells = which(deltas > qchisq(alpha, 1))
        #cat("lars badcells RSS dim:", badCells,"\n")
        badCells = which(lassOut$RSS > qchisq(alpha, F.rk))

        
        # check bad cells and impute 
        if(length(badCells) > 0){
            
            # corrupted cells indices
            ind     = min(max(badCells)+1,length(lassOut$nnzero.path))
            badinds = lassOut$nnzero.path[[ind]] # take at most the allowed one 
            #badinds = lassOut$ordering[seq_len(min(max(badCells)))]
            
            
            #  if(length(badinds) == p){
            #      cond.mhds = (x - mu) / sqrt(diag(Sigma))
            #  }else{
            # 
            # 
            #    # compute conditional covariance
            #    con.cov = ???==????
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
            # badinds = which(abs(sdres) > sqrt(qchisq(alpha, 1)))

          
            # impute data 
            if(length(badinds) == p){
                              
                Ximp[i, ]     = mu
                Zres.num[i, ] = (x - mu)
                Zres.den[i, ] = sqrt(diag(Sigma))
                Zres[i, ]     = Zres.num[i, ]/Zres.den[i, ]
            
            }else{
              

                # get best prediction
                best.guess = beta.fun(badinds, x, mu, precision)[badinds]
               
                
                # impute
                ximp          = X.new[i,]
                ximp[badinds] = ximp[badinds] - best.guess 
                Ximp[i, ]     = ximp
                
                # calculate residuals and Zres.num/Zres.den
                residual             = X.new[i, ] - ximp
                Zres.num[i, badinds] = residual[badinds]
                val                  = ginv(precision[badinds,badinds])    # inverse of subprecision matrix 
                Zres.den[i, badinds] = sqrt(diag(abs(val)))  
                Zres[i, badinds]     = Zres.num[i, badinds] / Zres.den[i,badinds]

            }
        }
    }
    
    # get outlying cells indices 
    indCells = which(abs(Zres)!=0)
    #indCells = which(abs(Zres) > sqrt(qchisq(alpha, 1)))

    # names
    colnames(Zres) = colnames(X); rownames(Zres) = rownames(X)
    colnames(Ximp) = colnames(X); rownames(Ximp) = rownames(X)

    # add back projected part
    Ximp.orig = Ximp + X.proj
        
    return(list(Ximp = Ximp, Ximp.orig = Ximp.orig, indCells = indCells, 
        Zres = Zres, Zres.den = Zres.den))
}


# function to compute possible paths of outlying cell blocks and 
# the corresponding mahalanobis distances 
dist.paths.fun <- function(X.new, grp.df, mu, predictors, precision, weights){
       
      # init
      n = dim(X.new)[1]
      p = dim(X.new)[2]
      
      # get guesses for which blocks to impute
      var = lapply(seq_len(n), function(idx){
       
      x            = X.new[idx, ]
      response     = predictors %*% (x - mu)
      weights.loc  = weights[idx,]
      
      # use groupwise lasso
      lassOut = lasso.fun(predictors, response, weights.loc, grp.df)
      lassOut$beta = beta.fun(lassOut$nnzero.path, x, mu, precision)
      lassOut$RSS  = RSS.fun(lassOut$beta, response, predictors)
      
      #cat("lass: ", lassOut$ordering, "\n")
      #cat("lass: ", round(lassOut$RSS,2),"\n")
      
      return(list(mhds = lassOut$RSS, 
                  nnzero.path = lassOut$nnzero.path,
                  added.vars  = lassOut$addedvars)) # mahalanobis distance 
      })
      
      # get maximum length and pad mahalanobis distances with NAS
      max.len = max(do.call("c",lapply(var, function(el){ length(el$mhds) })))
      mhds = lapply(var, function(el){ 
        el = el$mhds
        if(length(el) == max.len){ 
          return(el) 
        }else{  
          return(c(el, rep(NA, max.len - length(el))))
        } })
      mhds = do.call("rbind",mhds)
     
      # compute the differences in mhds
      mhds.diffs = t(apply(mhds, 1, function(el){ abs(diff(el)) }))
       
      # get nnzero paths
      nnzero.paths = lapply(var, function(el){
        return(el$nnzero.path)
      })
      
      # get added variables 
      added.vars =  lapply(var, function(el){
        return(el$added.vars)
      })
      
      # compute nbr of added variables in each step
      added.vars.len = lapply(added.vars, function(el){
        el = Reduce("c",lapply(el, function(subel){ length(subel) }))
        if(length(el) == dim(mhds.diffs)[2]){ 
          return(el) 
        }else{  
          return(c(el, rep(NA, dim(mhds.diffs)[2] - length(el))))
        } })
      added.vars.len = do.call("rbind", added.vars.len)
      
      return(list(mhds = mhds,
                  mhds.diffs = mhds.diffs,
                  nnzero.paths = nnzero.paths,
                  added.vars   = added.vars,
                  added.vars.len = added.vars.len))
}


# function that returns a matrix of flagged cells
flag.cells <- function(X.new,
                       mhds.diffs, 
                       nnzero.paths,
                       added.vars.len,
                       added.vars,
                       max.cells,
                       quant){
  
     # init
     n = dim(X.new)[1]
     p = dim(X.new)[2]
  
     #  take the biggest mahalanobis distances such that at most maxCol * n flagged 
     ties           = t(apply(mhds.diffs, 1, function(y)order(y))) + seq_len(n)*dim(mhds.diffs)[2]
     mhds.diffs.ord = order(mhds.diffs, ties, decreasing = TRUE)      # nas are ordered at the very end 
     #mhds.diffs.ord = unique(c(which(na.mask == 1), mhds.diffs))      # put nas in front  
     
     flagnbr.cols = rep(0, p)     # nbr of cells to impute for each column 
     drop.row     = rep(0, n)
     W            = matrix(0, n, p) # matrix of 1 and zeros
    
      # get outlying cells to impute 
      for(i in seq_len(length(mhds.diffs.ord))){
        
          idx        = mhds.diffs.ord[i]
          mhd.diff   = mhds.diffs[idx]        # current mahalanobis distance 
          row.idx    = (idx - 1) %% n + 1     # current row index in mhds matrix
          col.idx    = ((idx - 1) %/% n) + 1  # current column index in mhds matrix
          df.loc     = added.vars.len[idx]

          if(is.na(df.loc)){ next }
          threshold = qchisq(quant, df = df.loc)  # local threshold with degrees of freedom

          
          if(mhd.diff > threshold && !is.na(mhd.diff)){
            
            curr.path   = nnzero.paths[[row.idx]] # the current path of guesses
            path.idx    = min(col.idx + 1, length(curr.path))  # index in path of blocks
            curr.block  = curr.path[[path.idx]]
            added.block = added.vars[[row.idx]][[col.idx]]

            # check if total flagged cells would exceed allowed nbr 
            if(all((flagnbr.cols[added.block] + 1) <= max.cells) && drop.row[row.idx] == 0){
              
            #cat("row index: ",row.idx,"added block:", curr.block, "\n")  
                            
              W[row.idx, curr.block]     = 1
              flagnbr.cols[added.block]  = flagnbr.cols[added.block] + 1

            }else{ drop.row[row.idx] = 1 }
          }
      }
        
     
  return(W)
}


# function to calculate weights df
weights.fun <- function(X.new, mu, precision, grp.df){
  
  # init
  groups   = unique(grp.df$groups)
  n.groups = length(groups)
  
  # compute the precision matrices for each grouping 
  prec.groups = lapply(seq_len(n.groups), function(grp.idx){
     var.idx = grp.df$vars[grp.df$groups == groups[grp.idx]]
     precision[var.idx,var.idx]
  })
  
  # compute the weight estimates - for each group
  weights = sapply(seq_len(n.groups), function(grp.idx){
    var.idx = grp.df$vars[grp.df$groups == groups[grp.idx]]
    xs      = scale(X.new[,var.idx], center = mu[var.idx], scale = F)
    
    # mahalanobis distance
    mhd = sqrt(diag(xs %*% prec.groups[[grp.idx]] %*%  t(xs)))
    pmin(1, 1.5 / mhd) #1.5 / mhd #
  }, simplify = "matrix")
  
  return(weights)
}


singular.DI <- function(X, 
                        F.mat, 
                        grp.df, 
                        maxCol = 0.25,
                        alpha = 0.99, 
                        maxiter = 10,
                        conv.eps = 1e-6){
  
  
  #grp.df = data.frame(vars = c(1:p,4,5), groups = c(1,1,3,2,2,1,1))
  # grp.df = data.frame(vars = c(1:p), groups = 1:p)
  # 
  # maxCol = 0.25
  # alpha = 0.99
  # maxiter = 20
  # conv.eps = 1e-6

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
  F.svd  = F.mat.svd(F.mat)
  V.mat  = F.svd$V.mat
  V2.mat = V.mat %*% t(V.mat)
  F.rk   = dim(V.mat)[2]
  
  # orthogonal contribution
  X.proj = X %*% (diag(rep(1,p))-V2.mat)
  
  # replace NAs
  na.mask = is.na(X) + 0
  mu      = colMedians(X)
  for(i in seq_len(n)){
      indNA      = which(na.mask[i, ] == 1)
      X[i,indNA] = mu[indNA]
  }
  
  
  # calculate 
  # X.proj = X %*% (diag(rep(1,p))-V2.mat)
  # X.new  = X %*% V2.mat
  # 
  # # compute rescaling etc
  # res.obj   = rescaling(X.new, F.svd)
  # x.centers = res.obj$x.centers
  # X.new     = res.obj$x.mat.new
  # F.new     = F.mat#res.obj$F.mat.new
  # 
  # 
  # # compute initial mean/cov estimates
  # init.sol = cellWise:::DDCWcov(X.new)$Z
  # B.mat    = cov(init.sol %*% V.mat)             # check initial estimates again!!
  # mu       = colMeans(init.sol %*% V2.mat)       # check initial estimates again 
  # Sigma    = V.mat %*% B.mat %*% t(V.mat)        # check if best 
  # B.out    = red.ginv(B.mat)
  # precision = V.mat %*% B.out$Inv %*% t(V.mat)
  # predictors = V.mat %*% B.out$InvSqrt %*% t(V.mat)

  
  #######  or do 
  # compute rescaling etc
  res.obj   = rescaling(X, F.svd)
  x.centers = res.obj$x.centers
  X.new     = res.obj$x.mat.new
  F.new     = F.mat #res.obj$F.mat.new

  # compute initial mean/cov estimates
  init.sol = cellWise:::DDCWcov(X.new)$Z
  B.mat    = cov(init.sol %*% V.mat)             # check initial estimates again!!
  mu       = colMeans(init.sol %*% V2.mat)       # check initial estimates again 
  Sigma    = V.mat %*% B.mat %*% t(V.mat)        # check if best 
  B.out    = red.ginv(B.mat)
  precision = V.mat %*% B.out$Inv %*% t(V.mat)
  predictors = V.mat %*% B.out$InvSqrt %*% t(V.mat)
  
  # calculate 
  X.proj = X %*% (diag(rep(1,p))-V2.mat)
  X.new  = X %*% V2.mat
  #######

  # initialize vars for loop
  iter  = 0
  error = Inf
  max.cells = floor(maxCol * n)

  # run loop
  while(iter < maxiter && error > conv.eps){
    
     
     # compute current weights
     weights = weights.fun(X.new, mu, precision, grp.df)
     
     
     # calculate for each sample the mahalanobis distances 
     # of the given possible outlieying cell blocks 
     dps  = dist.paths.fun(X.new, grp.df, mu, predictors, precision, weights)
     mhds           = dps$mhds
     mhds.diffs     = dps$mhds.diffs
     nnzero.paths   = dps$nnzero.paths
     added.vars     = dps$added.vars
     added.vars.len = dps$added.vars.len
     
     # flag cells
     W.old = W

     W = flag.cells(X.new, mhds.diffs, nnzero.paths,
                       added.vars.len, added.vars,
                       max.cells, quant)
    print(norm(W-W.old,"O"))

    

     # impute data and calculate bias correction
     betas      = matrix(0, n, p)
     B.bias     = array(0,dim = dim(B.mat))
     B.bias.loc = array(0,dim = dim(B.mat))
     bias       = matrix(0, p, p)

     for(row.idx in seq_len(n)){
       
         imp.ind = which(!!W[row.idx,])
         p.imp   = length(imp.ind)
         if(length(imp.ind)!=0){
           
           x               = X.new[row.idx,]
           betas[row.idx,] = beta.fun(imp.ind, x, mu, precision)
           #B.bias.loc = Sigma[imp.ind,imp.ind] - Sigma[imp.ind,-imp.ind] %*% ginv(Sigma[-imp.ind,-imp.ind]) %*% Sigma[-imp.ind,imp.ind]#ginv(precision[imp.ind,imp.ind])
           B.bias.loc = ginv(precision[imp.ind,imp.ind])
           B.bias.loc = t(matrix(V.mat[imp.ind,], nrow = p.imp)) %*% B.bias.loc %*%  matrix(V.mat[imp.ind,], nrow = p.imp)
           B.bias     = B.bias + B.bias.loc
           
         }
     }
     
    # impute data 
    X.imp = X.new - betas
    #X.imp = X.imp %*% V2.mat
    
    # save current estimates
    mu.old    = mu
    B.mat.old = B.mat 
    Sigma.old = Sigma
    
    # update estimates 
    mu     = colMeans(X.imp %*% V2.mat)
    B.mat  = cov(X.imp %*% V.mat) 
    B.mat  = B.mat + 1 / (n-1) * B.bias                 # adding bias !
    Sigma  = V.mat %*% B.mat %*% t(V.mat)
    B.out  = red.ginv(B.mat)
    precision = V.mat %*% B.out$Inv %*% t(V.mat)
    predictors = V.mat %*% B.out$InvSqrt %*% t(V.mat)
  
    
    # update loop variables
    error = norm(B.mat - B.mat.old,"F") / (1 + norm(B.mat.old,"F"))
    error = max(error, norm(mu - mu.old,"2") / (1 + norm(mu.old,"2")))
    iter  = iter + 1
    
    cat("iter: ", iter, "error: ", error, "\n")
  }
  
  #  final mean and B.mat esitmates
  X.imp = X.imp %*% V2.mat
  X.imp = scale(X.imp, center = -x.centers, scale = F)   # add back mean
  mu.final    = colMeans(X.imp %*% V2.mat)
  B.mat.final = cov(X.imp %*% V.mat) 
  Sigma.final = V.mat %*% B.mat.final %*% t(V.mat)
  
  # final model estimates 
  est.final      = singular.cellImputation(X, V.mat, grp.df, mu.final, B.mat.final, alpha)
  corr.cells     = which(est.final$Zres != 0)
  corr.cells.arr = which(est.final$Zres != 0, arr.ind = T)
  corr.cells.arr = corr.cells.arr[order(corr.cells.arr[,1]),]
  colnames(corr.cells.arr) = c("n.idx", "p.idx")
  
  # go back 
  X.imp.proj = est.final$Ximp %*% V2.mat
  X.imp      = est.final$Ximp + X.proj
  Zres.proj  = est.final$Zres %*% V2.mat
  Zres.orig  = est.final$Zres

  
  return(list(X.imp          = X.imp,
              X.imp.proj     = X.imp.proj,
              corr.cells     = corr.cells, 
              corr.cells.arr = corr.cells.arr,
              Zres.proj      = Zres.proj,
              Zres.orig      = Zres.orig,
              mu             = mu.final,
              B.mat          = B.mat.final,
              Sigma          = Sigma.final))
}


# runs cross validation on given alphas

cv.alpha <- function(X, alphas = c(0.9,0.95,0.99)){
  
}






  x.corr.mat = log(robCompositions::expendituresEU)
  x.corr.mat = as.matrix(x.corr.mat)
  p = dim(x.corr.mat)[2]
 
  F.mat = diff.mat.fun(p,combn(1:p,2),rep(1,p*(p-1)/2))

  grp.df = cbind(1:12,1:12)
  #grp.df = cbind(1:5,1:5)
  grp.df = data.frame(grp.df)
  names(grp.df) = c("vars","groups")
  

  bla = cellHandlerExt(x.corr.mat, F.mat, grp.df, 0.99)
  plot(abs(!!bla$Zres))
  plot(abs(bla$Zres))
  
  library(reshape2)
  library(ggplot2)
  mod1 = cellWise::DI(x.corr.mat)
  bla = dim(melt(x.corr.mat))[1]
  plotdata = rbind(as.matrix(melt(x.corr.mat)),melt(mod1$Ximp))
  plotdata = cbind(plotdata, c(rep(1,bla),rep(2,bla)))
  colnames(plotdata) = c("Var1","Var2","value","imp")
  
  ggplot(plotdata,aes(x = Var2,y = value,col = as.factor(imp)))+ 
    geom_point() + facet_wrap(~Var1)+ylab("")+
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

  
  