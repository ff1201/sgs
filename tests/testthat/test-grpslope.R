#grpSLOPE and sgs do standardisations differently. sgs follows the convention of glmnet and SLOPE. Hence, no test for this.
test_that("solution reduces to grpslope when alpha=0 for uneven groups with no standardisation or intercept", {
  skip_if_not_installed("grpSLOPE")
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  n = nrow(X)
  lambda=1.5

  sgs = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda/n,   gFDR=0.1,standardise="none", intercept=FALSE) # lambda is scaled as grpSLOPE doesn't use 1/n factor 
  grpslope = grpSLOPE::grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda ="mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  grpslope_cost = sgs_convex_opt(X=X,y=y,beta=grpslope$beta,num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(grpslope_cost,sgs_cost,tol=1e-3)
})

test_that("solution reduces to grpslope when alpha=0 for even groups with no standardisation or intercept", {
  skip_if_not_installed("grpSLOPE")
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  n = nrow(X)
  lambda=1.5

  sgs = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda/n,   gFDR=0.1,intercept=FALSE,standardise="none") # lambda is scaled as grpSLOPE doesn't use 1/n factor
  grpslope = grpSLOPE::grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda = "mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  grpslope_cost = sgs_convex_opt(X=X,y=y,beta=grpslope$beta,num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta),
    tol = 1e-3
  )
  expect_equivalent(grpslope_cost,sgs_cost,tol=1e-3)
})

test_that("solution reduces to grpslope when alpha=0 for uneven groups with intercept but no standardisation", {
  skip_if_not_installed("grpSLOPE")
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  # grpslope implements the intercept differently, so centering y and X as this is done in sgs when intercept=TRUE
  X_center = apply(X,2,mean)
  X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
  y = data$y
  y=y-mean(y)
  n = nrow(X)
  lambda=1.5
  
  sgs = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda/n,   gFDR=0.1,intercept=TRUE,standardise="none") # lambda is scaled as grpSLOPE doesn't use 1/n factor
  grpslope = grpSLOPE::grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda ="mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta)[-1],
    tol = 1e-3
  )
})

test_that("solution reduces to grpslope when alpha=0 for even groups with intercept but no standardisation", { 
  skip_if_not_installed("grpSLOPE")
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  # grpslope implements the intercept differently, so centering y and X as this is done in sgs when intercept=TRUE
  X_center = apply(X,2,mean)
  X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
  y = data$y
  y=y-mean(y)
  n = nrow(X)
  lambda=1.5

  sgs = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda/n,   gFDR=0.1,intercept=TRUE,standardise="none") # lambda is scaled as grpSLOPE doesn't use 1/n factor
  grpslope = grpSLOPE::grpSLOPE(X=X, y=y, group = groups, fdr=0.1, lambda ="mean", max.iter=5000, normalize=FALSE,orthogonalize=FALSE,sigma=lambda)

  expect_equivalent(grpslope$beta,
    as.matrix(sgs$beta)[-1],
    tol = 1e-3
  )
})

test_that("sgs solution reduces to gslope when alpha=0 for uneven groups with no standardisation or intercept", {
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda=0.8

  gslope = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, standardise="none", intercept=FALSE)  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, alpha=0, standardise="none", intercept=FALSE)  
  
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  grpslope_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(gslope$beta),num_obs=1,gslope_seq=gslope$pen_gslope,slope_seq=gslope$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(grpslope_cost,sgs_cost,tol=1e-3)
})

test_that("sgs solution reduces to gslope when alpha=0 for even groups with no standardisation or intercept", {
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda=0.8

  gslope = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, intercept=FALSE,standardise="none") 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, alpha=0, standardise="none", intercept=FALSE)  
  
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=1,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  grpslope_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(gslope$beta),num_obs=1,gslope_seq=gslope$pen_gslope,slope_seq=gslope$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )
  expect_equivalent(grpslope_cost,sgs_cost,tol=1e-3)
})

test_that("sgs solution reduces to gslope when alpha=0 for uneven groups with intercept but no standardisation", {
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda=0.8
  
  gslope = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, intercept=TRUE,standardise="none") 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, alpha=0, standardise="none", intercept=TRUE)  
  
  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )
})

test_that("sgs solution reduces to gslope when alpha=0 for even groups with intercept but no standardisation", {
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda=0.8

  gslope = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, intercept=TRUE,standardise="none") 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, alpha=0, standardise="none", intercept=TRUE)  
  
  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )
})

test_that("sgs solution reduces to gslope when alpha=0 for uneven groups with intercept and standardisation", {
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda=0.8

  gslope = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, intercept=TRUE,standardise="l2") 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, alpha=0, standardise="l2", intercept=TRUE)  
  
  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )
})

test_that("sgs solution reduces to gslope when alpha=0 for even groups with intercept and standardisation", {
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda=0.8

  gslope = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, intercept=TRUE,standardise="l2") 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, alpha=0, standardise="l2", intercept=TRUE)  
  
  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )
})

test_that("sgs solution reduces to gslope when alpha=0 under logistic regression", {
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  n = 50
  p = 100
  X = as.matrix(rnorm_multi(n=n,vars=p,mu=0,sd=1,r=0))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  lambda=0.01

  gslope = fit_gslope(X=X,y=y, groups=groups, type="logistic", lambda=lambda, gFDR=0.1, standardise="l2", intercept=FALSE)  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="logistic", lambda=lambda, gFDR=0.1, alpha=0, standardise="l2", intercept=FALSE)  
  
  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

})

test_that("sgs solution reduces to gslope when alpha=0 under logistic regression with no screening", {
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  num_groups = length(unique(groups))

  n = 50
  p = 100
  X = as.matrix(rnorm_multi(n=n,vars=p,mu=0,sd=1,r=0))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  lambda=0.8

  gslope = fit_gslope(X=X,y=y, groups=groups, type="logistic", lambda=lambda, gFDR=0.1, standardise="l2", intercept=FALSE, screen = FALSE)  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="logistic", lambda=lambda, gFDR=0.1, alpha=0, standardise="l2", intercept=FALSE, screen = FALSE)  
  
  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

})

test_that("sgs solution reduces to gslope when alpha=0 under linear regression with no screening", {
  groups = rep(1:20,each=5)
  num_groups = length(unique(groups))

  data = gen_toy_data(p=100,n=50,rho = 0,seed_id = 3,grouped = TRUE,group_sparsity=0.1,groups = groups,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda=0.8

  gslope = fit_gslope(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, intercept=TRUE,standardise="l2", screen = FALSE) 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, gFDR=0.1, alpha=0, standardise="l2", intercept=TRUE, screen = FALSE)  
  
  expect_equivalent(as.matrix(gslope$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )
})