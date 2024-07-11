# gglasso doesn't apply standardisation, so can't compare directly with standardisation
test_that("solution reduces to grplasso when alpha=0 and constant weights, with uneven groups, with no intercept and standardisation", {
  skip_if_not_installed("gglasso")
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  p=100
  n=50
  data = gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1, alpha=0, vFDR=0.1,standardise="none", gFDR=0.1,intercept=FALSE,v_weights=rep(0,p),w_weights=rep(lambda,length(table(groups))))
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=FALSE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))
   
  gglasso_cost = sgs_convex_opt(X=X,y=y,beta=gglasso$beta,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(gglasso$beta,
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(gglasso_cost,sgs_cost,tol=1e-3)
})
  
test_that("solution reduces to grplasso when alpha=0 and constant weights, with even groups, with no intercept and standardisation", {
  skip_if_not_installed("gglasso")
  groups = rep(1:20,each=5)
  p=100
  n=50
  data = gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1, standardise="none", alpha=0, vFDR=0.1, gFDR=0.1,intercept=FALSE,v_weights=rep(0,p),w_weights=rep(lambda,length(table(groups))))
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=FALSE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))
   
  gglasso_cost = sgs_convex_opt(X=X,y=y,beta=gglasso$beta,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(gglasso$beta,
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(gglasso_cost,sgs_cost,tol=1e-3)
})

test_that("solution reduces to grplasso when alpha=0 and constant weights, with uneven groups, with intercept but no standardisation", {
  skip_if_not_installed("gglasso")
  groups = c(rep(1:5, each=3),
           rep(6:11, each=4),
           rep(12:16, each=5), 
           rep(17:22,each=6))
  p=100
  n=50
  data = gen_toy_data(p=p,n=n,rho = 0,seed_id = 4,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1, alpha=0, vFDR=0.1,standardise="none", gFDR=0.1,intercept=TRUE,v_weights=rep(0,p),w_weights=rep(lambda,length(table(groups))))
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=TRUE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))
   
  gglasso_cost = sgs_convex_opt(X=X,y=y,beta=c(gglasso$b0, gglasso$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(c(gglasso$b0, gglasso$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(gglasso_cost,sgs_cost,tol=1e-3)
})
  
test_that("solution reduces to grplasso when alpha=0 and constant weights, with even groups, with intercept but no standardisation", {
  skip_if_not_installed("gglasso")
  groups = rep(1:20,each=5)
  p=100
  n=50
  data = gen_toy_data(p=p,n=n,rho = 0,seed_id = 100,grouped = TRUE,groups = groups,group_sparsity=0.2,var_sparsity=1,orthogonal = FALSE)
  X = data$X
  y = data$y
  lambda = 0.8
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1,standardise="none", alpha=0, vFDR=0.1, gFDR=0.1,intercept=TRUE,v_weights=rep(0,p),w_weights=rep(lambda,length(table(groups))))
  gglasso = gglasso::gglasso(x=X,y=y, group=groups, lambda=lambda,intercept=TRUE,eps=1e-5, loss="ls",pf=sqrt(table(groups)))
   
  gglasso_cost = sgs_convex_opt(X=X,y=y,beta=c(gglasso$b0, gglasso$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(c(gglasso$b0, gglasso$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(gglasso_cost,sgs_cost,tol=1e-3)
})
  