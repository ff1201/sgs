### Comparison to seagull package: seagull uses a proximal algorithm. It implements intercept and standardisation together, so can't compare individually
test_that("solution reduces to sgl when using constant weights, with uneven groups, with no intercept and standardisation", {
  library(seagull)
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  alpha = 0.3
  sgl = seagull(y=y,X=matrix(0,nrow=dim(X)[1],ncol=dim(X)[2]),groups=groups,alpha=alpha,Z=X,max_lambda=lambda,loops_lambda=1,standardize = FALSE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  sgl_cost = sgs_convex_opt(X=X,y=y,beta=t(sgl$random_effects),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  
  expect_equivalent(t(sgl$random_effects),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(sgl_cost,sgs_cost, tol=1e-3)
})

test_that("solution reduces to sgl when using constant weights, with even groups, with no intercept and standardisation", {
  library(seagull)
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 5,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  alpha = 0.3
  sgl = seagull(y=y,X=matrix(0,nrow=dim(X)[1],ncol=dim(X)[2]),groups=groups,alpha=alpha,Z=data$X,max_lambda=lambda,loops_lambda=1,standardize = FALSE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  sgl_cost = sgs_convex_opt(X=X,y=y,beta=t(sgl$random_effects),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  
  expect_equivalent(t(sgl$random_effects),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(sgl_cost,sgs_cost, tol=1e-3)
})
 
### Comparison to SGL package: SGL uses a descent algorithm. SGL appears to standardise differently, so no comparison made.
test_that("solution reduces to sgl when using constant weights, with even groups, with no intercept and no standardisation", { 
  library(SGL)
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  # Can't turn off intercept for SGL, so removing it from data
  y = y - mean(y)
  X_center = apply(X,2,mean)
  X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
  lambda=0.8
  alpha = 0.3
  sgl = SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  sgl_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgl$beta[,2]),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  
  expect_equivalent(sgl$beta[,2],
    as.matrix(sgs$beta),
    tol = 1e-3
  )
  expect_equivalent(sgl_cost,sgs_cost, tol=1e-3)
})

test_that("solution reduces to sgl when using constant weights, with uneven groups, with no intercept and no standardisation", { 
  library(SGL)
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  # Can't turn off intercept for SGL, so removing it from data
  y = y - mean(y)
  X_center = apply(X,2,mean)
  X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
  lambda=0.8
  alpha = 0.3
  sgl = SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  sgl_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgl$beta[,2]),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=FALSE)
  
  expect_equivalent(sgl$beta[,2],
    as.matrix(sgs$beta),
    tol = 1e-3
  )
  expect_equivalent(sgl_cost,sgs_cost, tol=1e-3)
})
 
test_that("solution reduces to sgl when using constant weights, with even groups, with intercept but no standardisation", { 
  library(SGL)
  n = 50
  p = 100
  groups = rep(1:20,each=5)
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  alpha = 0.3
  sgl = SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(c(sgl$beta[,2]), # SGL seems to calculate the intercept different to other packages
    as.matrix(sgs$beta[-1]),
    tol = 1e-3
  )
})

test_that("solution reduces to sgl when using constant weights, with even groups, with intercept but no standardisation", { 
  library(SGL)
  n = 50
  p = 100
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 10,grouped = TRUE, groups=groups,group_sparsity=0.3,var_sparsity=0.5,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  alpha = 0.3
  sgl = SGL(list(x=X,y=y), index=groups, type = "linear",nlam=1,lambdas=c(0,lambda),alpha=alpha,standardize=FALSE)
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=alpha, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="none",w_weights = rep(1,length(table(groups))),v_weights = rep(1,p))
  
  expect_equivalent(c(sgl$beta[,2]), # SGL seems to calculate the intercept different to other packages
    as.matrix(sgs$beta[-1]),
    tol = 1e-3
  )
})
