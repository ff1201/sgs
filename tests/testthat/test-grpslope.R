# grpSLOPE and sgs do standardisations differently. sgs follows the convention of glmnet and SLOPE. Hence, no test for this.
test_that("solution reduces to grpslope when alpha=0 for uneven groups with no standardisation or intercept", {
  options(warn=-1) # to surpress warning about using grpSLOPE package for alpha = 0
  library(grpSLOPE)
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  data = generate_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  n = dim(X)[1]
  lambda=1.5
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda/n, alpha=0, vFDR=0.1, gFDR=0.1,standardise="none", intercept=FALSE) # lambda is scaled as grpSLOPE doesn't use 1/n factor
   
  grpslope = grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda ="mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  grpslope_cost = sgs_convex_opt(X=X,y=y,beta=grpslope$beta,num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(grpslope_cost,sgs_cost,tol=1e-3)
})

test_that("solution reduces to grpslope when alpha=0 for even groups with no standardisation or intercept", {
  options(warn=-1) # to surpress warning about using grpSLOPE package for alpha = 0
  library(grpSLOPE)
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = generate_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  n = dim(X)[1]
  lambda=1.5
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda/n, alpha=0, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none") # lambda is scaled as grpSLOPE doesn't use 1/n factor
   
  grpslope = grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda = "mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  grpslope_cost = sgs_convex_opt(X=X,y=y,beta=grpslope$beta,num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(grpslope_cost,sgs_cost,tol=1e-3)
})

test_that("solution reduces to grpslope when alpha=0 for uneven groups with intercept but no standardisation", {
  options(warn=-1) # to surpress warning about using grpSLOPE package for alpha = 0
  library(grpSLOPE)
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  data = generate_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  # grpslope implements the intercept differently, so centering y and X as this is done in sgs when intercept=TRUE
  X_center = apply(X,2,mean)
  X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
  y = data$y
  y=y-mean(y)
  y = data$y
  n = dim(X)[1]
  lambda=1.5
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda/n, alpha=0, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="none") # lambda is scaled as grpSLOPE doesn't use 1/n factor
  grpslope = grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda ="mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta)[-1],
    tol = 1e-3
  )

})

test_that("solution reduces to grpslope when alpha=0 for even groups with intercept but no standardisation", { # not sure why there is a slight difference here
  options(warn=-1) # to surpress warning about using grpSLOPE package for alpha = 0
  library(grpSLOPE)
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = generate_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  # grpslope implements the intercept differently, so centering y and X as this is done in sgs when intercept=TRUE
  X_center = apply(X,2,mean)
  X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
  y=data$y
  y=y-mean(y)
  y = data$y
  n = dim(X)[1]
  lambda=1.5
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda/n, alpha=0, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="none") # lambda is scaled as grpSLOPE doesn't use 1/n factor
  grpslope = grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda ="mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta)[-1],
    tol = 1e-3
  )

})

