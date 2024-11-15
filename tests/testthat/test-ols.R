test_that("unregularized gaussian models reduces to OLS, with no intercept", {
  set.seed(3)
  n=200
  p=10
  X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y <- X %*%rnorm(10,mean=0,sd=sqrt(10)) + rnorm(200,mean=0,sd=1)

  groups = 1:p
  lm_fit = lm(y ~ as.matrix(X) - 1)
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=0, alpha=1, vFDR=0.1, gFDR=0.1, intercept=FALSE, standardise="none")

  ols_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(as.numeric(coef(lm_fit))),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=c(0,groups),intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=c(0,groups),intercept=FALSE)

  expect_equivalent(as.matrix(coef(lm_fit)),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(ols_cost, sgs_cost, tol=1e-3)
})

test_that("unregularized gaussian models reduces to OLS, with intercept", {
  set.seed(3)
  n=200
  p=10
  X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y <- X %*%rnorm(10,mean=0,sd=sqrt(10)) + rnorm(200,mean=0,sd=1)
  groups = 1:p
  lm_fit = lm(y ~ as.matrix(X))
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=0, alpha=1, vFDR=0.1, gFDR=0.1, standardise="none", intercept=TRUE)

  ols_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(as.numeric(coef(lm_fit))),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=c(0,groups),intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=c(0,groups),intercept=TRUE)

  expect_equivalent(coef(lm_fit),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(ols_cost, sgs_cost, tol=1e-3)
})