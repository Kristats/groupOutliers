library(limSolve)
library(genlasso)
library(ICSNP)
library(MASS)
library(lars)
library(cellWise)
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


# function to get the coordinate map (clr but for the D.F.matmat matrix)
coord.map <- function(x.mat, F.mat){
  
  # init
  const.med = 1.4826
  tol = sqrt(.Machine$double.eps)
  N  = dim(x.mat)[1]
  px = dim(x.mat)[2]
  pD = dim(F.mat)[1]
  
  # compute the eigenvecs
  D.svd     = svd(x = F.mat, nv = min(dim(F.mat)), nu = 0)
  positive  = min(sum(D.svd$d >  max(tol * D.svd$d, 0)), min( dim(F.mat)))
  x.trans = t(D.svd$v[,1:positive] %*% t(D.svd$v[,1:positive]) %*% t(x.mat))
  
  # names
  colnames(x.trans) = colnames(x.mat)
  rownames(x.trans) = rownames(x.mat)
  
  return(x.trans)
}


# rescaling function
rescaling <- function(x.mat, F.mat){
  
  # init
  const.med = 1.4826
  tol = sqrt(.Machine$double.eps)
  N  = dim(x.mat)[1]
  px = dim(x.mat)[2]
  pD = dim(F.mat)[1]
  
  # compute the eigenvecs
  D.svd     = svd(x = F.mat, nv = min(dim(F.mat)), nu = 0)
  positive  = min(sum(D.svd$d >  max(tol * D.svd$d, 0)), min(dim(F.mat)))
  x.trans.svd = t(diag((D.svd$d[1:positive])) %*% t(D.svd$v[,1:positive]) %*% t(x.mat))
   
  # compute robust center 
  x.centers = spatial.median(x.trans.svd) # in transformed coordinates
  x.centers = D.svd$v[,1:positive] %*% diag(1 / (D.svd$d[1:positive])) %*% x.centers
  
  # center data
  x.mat.new = scale(x.mat, center = x.centers, scale = F) # sapply(1:px, function(p.idx){ x.mat[,p.idx] - x.centers[p.idx]})
  colnames(x.mat.new) = colnames(x.mat)
  rownames(x.mat.new) = rownames(x.mat)
  
  # now compute scales
  Dx.scales = sapply(1:pD, function(pd.idx){ const.med * median(abs(F.mat[pd.idx,] %*% t(x.mat.new)))  }) 
  F.mat.new = diag(1 / Dx.scales) %*% F.mat
  
  
  return(list(x.mat.new = x.mat.new, F.mat.new = F.mat.new, x.centers = x.centers, Dx.scales = Dx.scales))
}


# get the right svd vectors corresponding to non-zero singular values
V.mat.fun <- function(F.mat){
  
  # init
  tol     = sqrt(.Machine$double.eps)
  p1      = dim(F.mat)[1] 
  p2      = dim(F.mat)[2]

  # compute the eigenvecs
  D.svd     = svd(x = F.mat, nv = min(dim(F.mat)), nu = 0)
  positive  = min(sum(D.svd$d >  max(tol * D.svd$d, 0)), min(dim(F.mat)))
  Dv.mat    = D.svd$v[,1:positive]
  Dv.mat.o  = matrix(D.svd$v[,-c(1:positive)], ncol = length(D.svd$d) - positive)
  
  return(list(V.mat = Dv.mat, V.mat.o = Dv.mat.o))
}


# calculate the precision matrix from from x and V.mat 
P.mat.fun <- function(x.mat, col.centers, V.mat){
  
  x.centered  = scale(x.mat, center = col.centers, scale = F)
  x.trans     = t(t(V.mat) %*%  t(x.centered))
  x.trans.inv = red.ginv(1 / (dim(x.trans)[1]-1) * t(x.trans) %*% x.trans )
  P.mat.sq = V.mat %*% x.trans.inv$InvSqrt %*% t(V.mat)
  P.mat    = V.mat %*% x.trans.inv$Inv %*% t(V.mat)
  S.mat    = V.mat %*% cov(x.trans) %*% t(V.mat)
  

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


# not constrained lasso lars with weights  
lars.ext.unconstr <- function(predictors, response, constrs, V2.mat.proj, weights){
  
  # set some coeffs to zero that are too small
  eps.zero  = 1e-8
  eps.ols   = 1e-2
  
  # init 
  p = dim(predictors)[2]
  
  # solve with genlasso 
  mod = lars(x = predictors %*% diag(1 / weights), y  = response, type = "lar", normalize = F, intercept = F)
  #mods = sapply(1:p,function(idx){mean((lmRob(y ~.-1, data = data.frame(x = rbind(predictors %*% diag(1 / weights),diag(rep(1e-2,p)))[,idx], 
  #      y  = c(response,rep(0,p))))$T.residuals)^2) })
  mod$coefs = t(mod$beta %*% diag( weights))
  mod$coefs[abs(mod$coefs) < eps.zero] = 0
  
   # reduce to get coefficient path 
  nnzero.path    = apply(mod$coefs, 2, function(p.vec){which(p.vec!=0)})
  nnzero.path    = c(nnzero.path, list(1:p))                              # add all variables as well 
  ordering       = Reduce(function(w,v){ unique(c(w,v)) }, nnzero.path)
  
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
  #beta = beta %*% V2.mat.proj 
  
  #  get RSS
  RSS = apply(beta, 1, function(u){ sum((response - predictors %*% u)^2) })
  
  
  return(list(RSS = RSS, beta = beta, ordering = ordering, nnzero.path = nnzero.path))
}

#least squares
gen.ols <- function(x, y, idx, eps){
 
  # get reduced constraints
  loc.coeffs = lm.fit(x = rbind(matrix(x[,idx], ncol = length(idx)), sqrt(eps) * diag(rep(1,length(idx)))), 
                          y = c(y, rep(0,length(idx))))$coefficients

  
  return(loc.coeffs) 
}



# function to calculate weights beforehand for all samples...
weights.fun <- function(x.mat, x.imp, F.mat){
  
  # init
  N     = dim(x.mat)[1]
  N.imp = dim(x.imp)[1]
  p     = dim(x.mat)[2]
  pD    = dim(F.mat)[2]
  
  # centers and scales for weigths of lars
  #x.Dmat.centers = sapply(1:pD, function(pD.idx){ spatial.median(t(F.mat[,pD.idx] %*% t(x.mat[,pD.idx]))) })
  x.Dmat.centers = sapply(1:pD, function(pD.idx){ rowMeans(F.mat[,pD.idx] %*% t(x.imp[,pD.idx])) })

  x.Dmat.scales  = sapply(1:pD, function(pD.idx){ 
    
    var.loc = sapply(1:N.imp, function(n.idx){ F.mat[,pD.idx] %*% t(x.imp[n.idx,pD.idx]) - x.Dmat.centers[,pD.idx] })
    var.loc =  1.48 * ( median(apply(var.loc, 2, function(u){ norm(u,"2") }))) # sqrt( 1/(N.imp-1) * sum(apply(var.loc, 2, function(u){ norm(u,"2")^2 }))) #
    var.loc
    
  })
  
  # compute outliers viar Lars
  weights.mat = lapply(1:N, function(n.idx){
    
  # get sample
  u.vec = x.mat[n.idx,]
    
  # calculate weights
  weights = sapply(1:pD, function(pD.idx){  
    w.loc = norm(as.numeric(u.vec[pD.idx]) * F.mat[,pD.idx] - x.Dmat.centers[,pD.idx],"2") 
    w.loc = w.loc / x.Dmat.scales[pD.idx]
  })
  weights = pmin(1, 1.5 / weights) # 
  weights = weights / sum(weights)
  
  })
  weights.mat = Reduce("rbind",weights.mat)
  #weights.mat = diag(1/rowSums(weights.mat)) %*% weights.mat
  
  return(weights.mat)  
}

# implementation of the extension of the cellHandler algorithm 
cellHandlerExt <- function(X, mu, F.mat, P.out, V2.mat.proj, constrs, quant = 0.99){
  
    # init
    eps.ols = 1e-2
    eps     = 1e-6
    n = dim(X)[1]
    d = dim(X)[2]

    
    # replace X by projected part 
    X.proj = X %*% constrs %*% t(constrs)
    X      = X %*% V2.mat.proj
    
    # get precision matrix / covariance matrix & predictors
    P.mat      = P.out$P.mat
    Sigma      = P.out$S.mat
    predictors = P.out$P.mat.sq
    scales     = sqrt(diag(Sigma))
    d.ker      = dim(constrs)[2]   # dimension of kernel 
        
    
    # compute weight matrix for this setting 
    #weights.mat = weights.fun(X, x.imp, F.mat)
    weights.mat = t(apply(X, 1, function(u){pmin(1, 1.5/ (abs(u -  mu) / scales)) }))
    
    # init next 
    Ximp     = X
    naMask   = is.na(X) + 0   ###### NOT USED YET IN LAROUT: TODO 
    Zres     = matrix(0, n, d)
    Zres.den = matrix(0, n, d)
    Zres.num = Zres.den
    for (i in seq_len(n)) {
        x <- X[i, ]
        indNA <- which(naMask[i, ] == 1)
        x[indNA] <- mu[indNA]
        response <- P.out$P.mat.sq %*% (x - mu)
        weights  <- weights.mat[i,]
        
        # maybe use these weights - still 
        #weights = cellWise:::huberweights(x = (x -  mu) / scales, b = 1.5)
        
        # use lars 
        #larOut = lars.ext(predictors, response, constrs, V2.mat.proj, lam.max, weights)
        larOut = lars.ext.unconstr(predictors, response, constrs, V2.mat.proj, weights)
       
        # path index of possible outlying cells 
        deltas   =  abs((diff(larOut$RSS)))
        badCells = which(deltas > qchisq(quant, 1))
       
        # extend if NAs in data
        if (length(indNA) > 0) { #### to look at 
            badCells = unique(indNA, badCells)
        }
        
        # check bad cells and impute 
        if (length(badCells) > 0) {
            
            # corrupted cells indices
            badinds <- larOut$ordering[seq_len(min(max(badCells),length(larOut$ordering)))]

            # impute corrupted data
            if(length(badinds) == d){
                sdres = (x - mu) / sqrt(diag(Sigma))
            }else if(length(badinds) == 1){
                sdres = rep(0, d)
                #stdsc = Sigma[badinds,badinds]
                #stdsc = stdsc - Sigma[badinds,-badinds]%*%  ginv(Sigma[-badinds,-badinds]) %*% Sigma[-badinds,badinds]
                #stdsc = sqrt(diag(abs(stdsc)))
                stdsc = sqrt(abs(P.mat[badinds,badinds])) 
                #  P.mat[badinds,badinds] %*% (X[i,badinds]- mu[badinds] - best.guess)
                #  -P.mat[badinds,-badinds] %*% (X[i,-badinds]- mu[-badinds]) 
                
                if(stdsc < 1e-4){ 
                  sdres[badinds] = 0 
                }else{
                  sdres[badinds] = abs(larOut$beta[max(badCells) + 1, badinds])
                  sdres[badinds] = sdres[badinds] / stdsc 
                }
                #print(paste("sample nbr:", i, "...denom:", round(stdsc, 4)))
                
            }else{
                sdres = rep(0, d)
                #stdsc = Sigma[badinds,badinds] - Sigma[badinds,-badinds]%*% ginv(Sigma[-badinds,-badinds])%*%Sigma[-badinds,badinds]
                stdsc =  sqrt(diag(abs(ginv(P.mat[badinds,badinds])))) # sqrt(diag(abs(stdsc))) 
                sdres[badinds] = abs(larOut$beta[max(badCells) + 1, badinds])
                sdres[badinds] = sdres[badinds] / stdsc
                
                #print(paste0("sample nbr:", i, "...denom:"))
                #print(round(stdsc, 4))
            }
            
            # check which really to replace after correcting for variance 
            badinds <- which(abs(sdres) > sqrt(qchisq(quant,1)))              # without this actually better !!
            
            if (length(indNA) > 0) {
                badinds <- unique(indNA, badinds)
            }
            
            
            # check if any left
            if(length(badinds) == 0){ next }
            
            # impute data 
            if (length(badinds) == d) {
                              
                Ximp[i, ]     = mu
                Zres.num[i, ] = (X[i, ] - mu)
                Zres.den[i, ] = sqrt(diag(Sigma))
                Zres[i, ]     = Zres.num[i, ]/Zres.den[i, ]
            
            }else{
              
                # fit constraint ols to get best descreasing subvector for MD
                # is the best guess for centered data 
                #best.guess = constr.ols(predictors, response, t(constrs), badinds)
                best.guess = gen.ols(predictors, response, badinds, eps.ols)
                
                # impute
                ximp          = X[i, ]
                ximp[badinds] = ximp[badinds] - best.guess 
                Ximp[i, ]     = ximp
                
                # calculate residuals and Zres.num/Zres.den
                residual             = X[i, ] - ximp
                Zres.num[i, badinds] = residual[badinds]
                #sdres = Sigma[badinds,badinds] - Sigma[badinds,-badinds]%*% ginv(Sigma[-badinds,-badinds])%*%Sigma[-badinds,badinds]
                #Zres.den[i, badinds] = sqrt(diag(abs(sdres)))  
                Zres.den[i, badinds] = sqrt(diag(abs(ginv(P.mat[badinds,badinds]))))
                Zres[i, badinds]     = Zres.num[i, badinds] / Zres.den[i,badinds]
                #print(paste0("sample nbr:", i, "...denom:"))
                #print(round( Zres.den[i, badinds], 4))
            }
        }
    }
    
    
    indCells = which(abs(Zres) > sqrt(qchisq(quant, 1)))
    indNAs   = which(naMask == 1)
    indCells = setdiff(indCells, indNAs)

    #names
    colnames(Zres) = colnames(X); rownames(Zres) = rownames(X)
    colnames(Ximp) = colnames(X); rownames(Ximp) = rownames(X)

    # add bac projected part
    Ximp.orig = Ximp + X.proj
        
    return(list(Ximp = Ximp, Ximp.orig = Ximp.orig, indCells = indCells, indNAs = indNAs, 
        Zres = Zres, Zres.den = Zres.den))
}


# the whole generalized DI algorithm 
cell.est <- function(x.mat, F.mat, maxCol = 0.25, quant = 0.99, conv.eps = 1e-5, maxits = 10){
  
  
  # TODO: bias correction for the degenerate case ??!
  
  # TODO:: REDGINV computes full svd !! -> change this to reduced one --> implement
  
  
  # implementation of generalized DI algorithm 
  # init
  n = dim(x.mat)[1]
  d = dim(x.mat)[2]
  
  # calculate V.mat and orthogonal 
  V.mat.full  = V.mat.fun(F.mat)
  V.mat       = V.mat.full$V.mat
  V2.mat.proj = V.mat %*% t(V.mat)
  constrs     = V.mat.full$V.mat.o  
  
  
  # replace data by projected data
  x.mat.orth = x.mat %*% (constrs)
  x.mat      = x.mat %*% V2.mat
  
  
  # center data and rescale F.mat
  res        = rescaling(x.mat, F.mat) # center x.mat and rescale columns of F.mat - doenst change kernel 
  x.centers  = res$x.centers
  x.mat.new  = res$x.mat.new
  F.mat.new  = res$F.mat.new
  
  # run DDC to get initial estimates
  #start.sol = DDC(X = as.matrix(x.mat.new %*% t(F.mat.new) %*% F.mat.new))$Ximp 
  #start.sol = DDC(X = x.mat.new)$Ximp  # DI(x.mat.new)# DI(as.matrix(x.mat.new %*% t(F.mat) %*% F.mat))$Ximp  #
  #start.sol = cellHandlerExt(x.mat.new, 
  #                             colMeans(start.sol), 
  #                             F.mat,
  #                             P.mat.fun(start.sol,colMeans(start.sol),V.mat),
  #                             V2.mat,V.mat.o,quant = quant)$Ximp
  start.sol = cellWise:::DDCWcov(x.mat.new)$Z
  
  # compute inital estimates
  mu     = rescaling(start.sol, F.mat.new)$x.centers #x.centers
  P.out  = P.mat.fun(start.sol, mu, V.mat)
  P.mat  = P.out$P.mat
  predictors = P.out$P.mat.sq
  Sigma      = P.out$S.mat
  scales     =  diag(P.out$S.mat)
  d.ker      = dim(constrs)[2]   # dimension of kernel 
        
  
  
  #x.mat.new = start.sol
  #n         = dim(x.mat.new)[1]
  

  # intialize for loop
  iter      =  0
  error     = Inf
  M         = floor(maxCol * n)
  x.imp.mat = x.mat.new
  naMask    = is.na(x.mat) + 0
  
  # replace na's
  for(i in seq_len(n)){
      indNA    = which(naMask[i, ] == 1)
      x.mat.new[i,indNA] = mu[indNA]
  }
  
  # all.weights = weights.fun(x.mat.new, F.mat.new)

  
  while((iter < maxits) && (error > conv.eps)){
  
        predictors = P.out$P.mat.sq
        betamat    = array(0, c(n, d + 1, d))     # matrix of coefficients of path (2nd & 3rd dim) for each sample (first dim)
        Bmat       = array(0, c(n, d + 1, d, d))
        orderings  = matrix(0, n, d)
        distances  = matrix(0, n, d + 1)
        deltas     = matrix(0, n, d)
        
        # compute the weights
        #all.weights = weights.fun(x.mat.new, , F.mat)

        mhds = matrix(0,n,p)
        # solve for each sample the according lars problem 
        for (i in seq_len(n)) {
          
            # get x, assign weights and response
            x        = x.mat.new[i, ]
            response = predictors %*% (x - mu)
            #weights  = all.weights[i,]
            weights = cellWise:::huberweights(x = (x - mu) / scales, b = 1.5)
        
            # solve with lars 
            larOut = lars.ext.unconstr(predictors, response, constrs, V2.mat.proj, weights)
            
            # set deltas 
            deltas.loc = abs(diff(larOut$RSS))
            deltas.loc = rev(cummax(rev(deltas.loc)))
            deltas[i, larOut$ordering] = deltas.loc
            mhds.loc = apply(larOut$beta, 1, function(u){ t(x-mu-u)%*%P.out$P.mat %*%(x-mu-u) })[-c(p+1)]
            mhds[i,larOut$ordering] = mhds.loc[-(p+1)]
            
            betamat[i, , ] = larOut$beta
            #Bmat[i, , , ]  = larOut$biasMat TODO!!
            distances[i, ] = larOut$RSS
            orderings[i, ] = larOut$ordering
        }
        
        print(sum(apply(x.mat.new, 1, function(u){ t(u-mu)%*%P.out$P.mat %*%(u-mu) })))
        
        #sum(mhds > qchisq(quant,p-d.ker))
        #plot(mhds>qchisq(quant,p-d.ker))
        plot(A)
        #length(intersect(which(mhds>qchisq(quant,p-d.ker)), ind.corr))
        plot(deltas)
        
        
        # reorder and only take the biggest deltas such that at most maxCol * n flagged 
        tiebraker    = t(apply(orderings, 1, function(y)order(y))) + seq_len(n)*d
        deltas.order = order(deltas, tiebraker, decreasing = TRUE)   # ties in deltas are broken according to tiebraker 
        deltas.order = unique(c(which(naMask == 1), deltas.order))
        nbr.imp.cols = rep(0, d)                                     # nbr of cells to impute for each column 
        nbr.flag.cols= rep(0, n)                                     # nbr of flagged cells of each row 
        droppedPaths = rep(0, n)                                     # 0 or 1 (1 if row has been locked)
        W = matrix(0, n, d)
        
        # actual flagging - go down the list of ordered deltas 
        for(i in seq_len(length(deltas.order))){
            idx   = deltas.order[i]; delta = deltas[idx]  # delta info 
            n.ind = (idx - 1)%%n + 1   # row index of current cell
            
            # check if delta is bigger than threshold and could be a outlying cell
            if(delta > qchisq(quant, 1)){ 
                if(!droppedPaths[n.ind]){
                  d.ind <- ((idx - 1)%/%n) + 1  # column index of current cell 
                  # check if current column has too many flagged cells already
                  if(nbr.imp.cols[d.ind] < M){
                    nbr.flag.cols[n.ind] = nbr.flag.cols[n.ind] + 1
                    nbr.imp.cols[d.ind] = nbr.imp.cols[d.ind] + 1
                    W[n.ind, d.ind]     = 1    # set according cell mask to 1 
                  }else{ droppedPaths[n.ind] = 1 }
                }
            }else{ droppedPaths[n.ind] = 1 }
        }
        
        # get best guesses for imputation 
        finalBetas = matrix(0, n, d)
        finalBias  = matrix(0, d, d)
        for (i in seq_len(n)) {
            finalBetas[i, ] = betamat[i, nbr.flag.cols[i]+1, ]
            inds            = !!W[i,]
            p.imp           = sum(inds)
            if(p.imp!=0){
              bias.loc = as.matrix(ginv(P.out$P.mat[inds,inds]))
              bias.loc = t(matrix(V.mat[inds,], nrow = p.imp)) %*% bias.loc %*%  matrix(V.mat[inds,], nrow = p.imp)
              bias.loc = V.mat %*% bias.loc %*% t(V.mat)
              finalBias = finalBias +bias.loc
            }
        }
        
        # impute 
        x.imp.mat = x.mat.new - finalBetas

        # save old mu and Sigma 
        mu.old    = mu
        Sigma.old = Sigma
        
        # new covariance, precision and sqrt precision matrix
        mu             = colMeans(x.imp.mat %*% V2.mat)
        Sigma          = cov(x.imp.mat %*% V2.mat) + 1 / (n-1) * finalBias
        S.sqinv.mat    = red.ginv(Sigma)       # inverse and sqrt of inverse 
        P.out$P.mat    = S.sqinv.mat$Inv
        P.out$P.mat.sq = S.sqinv.mat$InvSqrt
        scales         = diag(Sigma)

        # update error and nbr of iterations 
        error = norm(Sigma.old - Sigma, "F") / (1 + norm(Sigma.old, "F") )
        error = max(error, norm(mu.old - mu,"2") / norm(mu.old,"2"))  
        print(error)
        
        # check PRF scores
        print(PRFscores(which(finalBetas != 0), ind.corr))
        print(which(finalBetas != 0))
   
        iter = iter + 1      
   }   
  
    # add center back to data  
    x.imp.mat = scale(x.imp.mat, center = -x.centers, scale = F)
    
    # final center and covariance/precision estimate 
    mu.final    = colMeans(x.imp.mat)
    P.out.final = P.mat.fun(x.imp.mat, mu.final, V.mat)
    
    # compute last cellHandler esitmate 
    est.final = cellHandlerExt(x.mat, mu.final, F.mat.new, P.out.final, V2.mat.proj , constrs, quant = quant)
    
    # add projection on null space as well  
    x.imp.final.proj = est.final$Ximp
    x.imp.final.orig = est.final$Ximp + x.mat.orth %*% t(constrs)
    x.mat.orig       = x.mat + x.mat.orth %*% t(constrs)
    Zres.orig        = est.final$Zres
    Zres.proj        = est.final$Zres
    #indCells         = cbind(est.final$indcells %% n, est.final$indcells %/% n + 1)
    #indCells         = which(abs(Zres.orig) != 0, arr.ind = T)
    #colnames(indCells) = c("n.idx","p.idx")
    indCells        = which(Zres.orig!=0,)
    
    #return
    return(list(x.mat.orig = x.mat.orig,
                x.mat.proj = x.mat,
                F.mat = F.mat,
                quant = quant, 
                x.imp.orig = x.imp.final.orig,
                x.imp.proj = x.imp.final.proj,
                Zres.orig = Zres.orig,
                Zres.proj = Zres.proj,
                Zres.den = est.final$Zres.den,
                indCells = indCells))
}




library(reshape2)
library(ggplot2)
N = 30; p = 5
# create D matrix
F.mat = diff.mat.fun(p,combn(1:p,2),rep(1,p*(p-1)/2))

n.rep  = 200
gamma  = 2
gammas = c(1,2,3,4,5,6,7,8,9,10)

# quant = 0.90 seems to be best ?!?!?? -> maybe plot for different quants.? or vary quant between 0-1???
quant = 0.90
all.scores = matrix(0, ncol = 8, nrow = length(gammas)) 
iter = 1
for(gamma in gammas){
  scores = lapply(1:n.rep,function(idx){

      print(idx)
        
      outl.gen   = outl.sample(p,N,0.2, gamma, F.mat)
      x.corr.mat = outl.gen$x.corr.mat
      x.mat      = outl.gen$x.mat
      ind.corr   = outl.gen$ind.corr
      
      
      #x.corr.mat = x.corr.mat %*% ( V.mats$V.mat) %*%  t(V.mats$V.mat)
      V.mats  = V.mat.fun(F.mat)
      V.mat   = V.mats$V.mat
      V.mat.o = V.mats$V.mat.o
      V2.mat  = V.mat%*%t(V.mat)
      #x.imp   = cellWise::DDC(as.matrix(x.corr.mat %*% (V.mat) %*% t(V.mat)))$Ximp
      #x.corr = x.corr.mat %*% (V.mat) %*% t(V.mat)   #----> makes it even worse!!!!
      
      
      constrs = V.mat.o
      
      x.corr.mat = x.corr.mat + rt(N,30) %*% t(constrs) 
      #x.corr.mat = x.corr.mat %*% V2.mat
      
      
      #res        = rescaling(x.corr.mat  %*% V2.mat, F.mat) # center x.mat and rescale columns of F.mat - doenst change kernel 
      #x.centers  = res$x.centers
      #x.mat.new  = res$x.mat.new
      #F.mat.new  = res$F.mat.new
      x.mat.new = x.corr.mat
      F.mat.new = F.mat     
      
      x.imp = cellWise:::DDCWcov(as.matrix(x.mat.new %*% V2.mat))$Z
      mu    = apply(x.imp, 2, function(u){mean(u)})

      mod = cellHandlerExt(x.corr.mat, mu, F.mat.new, P.mat.fun(x.imp, mu, V.mat), V2.mat, V.mat.o, quant = quant)

      
      #apply(x.corr.mat,1,function(u){t(u-colMeans(x.imp)) %*%P.mat.fun(x.imp, mu, V.mat)$P.mat  %*% (u-colMeans(x.imp))  }) > qchisq(p = 0.99,df = p-1)
      
      x.mat.new.imp = mod$Ximp %*% V2.mat
      P.out = P.mat.fun(x.mat.new.imp, colMeans(x.mat.new.imp), V.mat)

       #  
       # A =P.out$S.mat
       # A.old = matrix(0,p,p)
       # x.im.old = x.mat.new.imp
       # for(i in 1:50){
       #    
       #     
       #   #res        = rescaling(mod$Ximp  %*% V2.mat, F.mat) # center x.mat and rescale columns of F.mat - doenst change kernel 
       #   #x.centers  = res$x.centers
       #   #x.mat.new.imp  = res$x.mat.new
       # 
       #    
       #  mod = cellHandlerExt(x.corr.mat%*%V2.mat, 
       #                         mu,
       #                         F.mat.new, 
       #                         P.out, 
       #                         V2.mat, 
       #                         V.mat.o, 
       #                         quant = quant)
       #    
       #  W = matrix(0,N,p)  
       #  W[mod$indCells] = 1
       #    
       #  finalBias  = matrix(0, p, p)
       #  for (i in seq_len(N)) {
       #      inds            = !!W[i,]
       #      p.imp           = sum(inds)
       #      if(p.imp!=0){
       #        bias.loc = as.matrix(ginv(P.out$P.mat[inds,inds]))
       #        bias.loc = t(matrix(V.mat[inds,], nrow = p.imp)) %*% bias.loc %*%  matrix(V.mat[inds,], nrow = p.imp)
       #        bias.loc = V.mat %*% bias.loc %*% t(V.mat)
       #        finalBias = finalBias +bias.loc
       #      }
       #  }  
       #  
       #   # new covariance, precision and sqrt precision matrix
       #  x.mat.new.imp  = mod$Ximp
       #  mu             = colMeans(mod$Ximp %*% V2.mat)
       #  Sigma          = cov(mod$Ximp %*% V2.mat) + 1 / (N-1) * finalBias
       #  S.sqinv.mat    = red.ginv(Sigma)       # inverse and sqrt of inverse 
       #  P.out$S.mat    = Sigma
       #  P.out$P.mat    = S.sqinv.mat$Inv
       #  P.out$P.mat.sq = S.sqinv.mat$InvSqrt
       #  scales         = diag(Sigma)  
       #    
       #  print(norm((cov(x.mat.new.imp)-A.old),"F")/(1+norm(A.old,"F")))
       #  A.old = cov(x.mat.new.imp)
       #  #print(mean((mod$Ximp - x.im.old)^2)/(1+mean((x.im.old)^2)))
       #  x.im.old = x.mat.new.imp
       #  print(mod$indCells)
       #  print(PRFscores(mod$indCells, ind.corr))
       #  }

      
      #mod = cell.est(x.corr.mat, F.mat.new, maxit = 20, quant = quant, maxCol = 0.25); mod$Ximp = mod$x.imp.orig
      #mod2 = DI(x.corr.mat, quant = quant) # 
      mod2 = cellWise::cellHandler(x.corr.mat, mu, cov(x.imp) + diag(rep(1e-2,p)) ) #

      ext.cellH.imp = as.matrix(mod$Ximp %*% V2.mat)
      cellH.imp     = as.matrix(mod2$Ximp %*% V2.mat)
      err.ext.cellH = norm(x.mat %*% V2.mat - ext.cellH.imp ,"F") / sqrt(N)
      err.cellH     = norm(x.mat %*% V2.mat - cellH.imp ,"F") / sqrt(N)
      
      #sort(mod$indCells)
      #mag.order  = ind.corr[order(abs(x.mat[ind.corr]-x.corr.mat[ind.corr]), decreasing = T)]
      #(x.mat-x.corr.mat)[mag.order]
      
      return(c(ext.cellH = PRFscores(mod$indCells, ind.corr), 
               cellH     = PRFscores(mod2$indcells, ind.corr),
               err.ext.cellH = err.ext.cellH ,
               err.cellH     = err.cellH ))
        
})
scores = Reduce("rbind", scores)
#scores
colMeans(scores)[-c(7:8)]
ggplot(data.frame(melt(scores[,-c(7:8)])),aes(x = Var2, y = value)) + geom_boxplot()
ggplot(data.frame(melt(scores[,c(7:8)])),aes(x = Var2, y = value)) + geom_boxplot()

all.scores[iter,] = colMeans(scores)
iter = iter + 1


}

colnames(all.scores) = names(colMeans(scores))
#saveRDS(all.scores, "all.scores.rds")
all.scores = readRDS("all.scores.rds")


ggplot(melt(all.scores[,1:6]), aes(x = Var1,y = value)) + geom_line() + facet_wrap(~Var2) + ylim(0,1)
ggplot(melt(all.scores[,7:8]), aes(x = Var1,y = value)) + geom_line() + facet_wrap(~Var2)



all.scores = readRDS("all.scores2.rds")
ggplot(melt(all.scores[,1:6]), aes(x = Var1,y = value)) + geom_line() + facet_wrap(~Var2) + ylim(0,1)
ggplot(melt(all.scores[,7:8]), aes(x = Var1,y = value)) + geom_line() + facet_wrap(~Var2)




all.scores = readRDS("all.scores3.rds")
ggplot(melt(all.scores[,1:6]), aes(x = Var1,y = value)) + geom_line() + facet_wrap(~Var2) + ylim(0,1)
ggplot(melt(all.scores[,7:8]), aes(x = Var1,y = value)) + geom_line() + facet_wrap(~Var2)




# run once 
outl.gen   = outl.sample(p,N,0.2,2,F.mat)
x.corr.mat = outl.gen$x.corr.mat
x.mat      = outl.gen$x.mat
ind.corr   = outl.gen$ind.corr

mod = cell.est(x.corr.mat, F.mat, maxit = 30)
mod2 = DI(x.corr.mat)

PRFscores(mod$indCells, ind.corr)
PRFscores(mod2$indcells, ind.corr)


library(pheatmap)
pheatmap(abs(x.corr.mat-x.mat) , cluster_cols = F, cluster_rows  = F)
pheatmap(abs(mod$Zres.orig) , cluster_cols = F, cluster_rows  = F)
pheatmap(abs(mod2$Zres) , cluster_cols = F, cluster_rows  = F)



V.mats = V.mat.fun(F.mat)
mod = cellHandlerExt(x.corr.mat,mean(x.mat),F.mat,P.mat.fun(x.mat,colMeans(x.mat),V.mats$V.mat), V.mats$V.mat%*%t(V.mats$V.mat),V.mats$V.mat.o )

PRFscores(mod$indcells, ind.corr)
PRFscores(mod2$indcells, ind.corr)



library(reshape2)
library(ggplot2)
ppdf1 = melt(t(coord.map(mod$x.mat.proj,F.mat)))
ppdf2 = melt(t(coord.map(mod$x.imp.proj,F.mat)))
ppdf1 = cbind(as.data.frame(ppdf1),rep("x.orig",dim(ppdf1)[1]))
ppdf2 = cbind(as.data.frame(ppdf2),rep("x.imputed",dim(ppdf2)[1]))

ppdf = data.frame(Var1 = c(ppdf2$Var1,ppdf1$Var1),
           Var2 = c(ppdf2$Var2,ppdf1$Var2),
           value = c(ppdf2$value,ppdf1$value),
           col = c(factor(rep("x.imputed",dim(ppdf1)[1])),factor(rep("x.orig",dim(ppdf1)[1]))))
ppdf = as.data.frame(ppdf)
#colnames(ppdf) = c("Var1","Var2","value","col")

p1  = ggplot(ppdf, aes(y = value,x = Var1,col = (col)), alpha = 0.4) + 
      facet_wrap(~Var2) + 
      geom_point() + 
      scale_colour_manual(values = c("x.orig"="black","x.imputed"="red"))+
      xlab("Coordinate")+
      ylab("Clr-values")

p1






# try with some data 
x = robCompositions::phd[,7:11]; x[x==0] = 0.5 * min(x[x!=0]);  rownames(x) = robCompositions::phd[,2]
x.mat = as.matrix(log(x))
p = dim(x)[2]
N = dim(x)[1]


# first try !!!!
x.mat[1,] = x.mat[1,] + rep(2,p)

# second try
x.mat[1,2] = x.mat[1,2] * 3

# second simply try clr
x.mat = x.mat %*% V2.mat.proj


















    
  # degrees of freedom
  dofs = lapply(nnzero.path, function(idx){ 

     
    if(length(idx)==0){ 
      return(p) 
    }else if( length(idx)== p){  
      return(0)
    }else{
      dof = p - sum(diag(predictors[,idx] %*% ginv(P.mat[idx,idx]) %*% t( predictors[,idx]) )) 
      return(dof)
    }
    })
  dofs = Reduce("c", dofs)
  print(dofs)
        
  
