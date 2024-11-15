test_that("unregularized binomial models reduces to glm without intercept", {
  set.seed(3)
  n=200
  p=10
  X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y <- 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y <- ifelse(y>0.5,1,0)

  groups = 1:p
  glm_fit = glm(y ~ as.matrix(X)-1,family="binomial")
  sgs = fit_sgs(X=X,y=y, groups=groups, type="logistic", lambda=0, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none")

  ols_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(as.numeric(coef(glm_fit))),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=c(0,groups), intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=c(0,groups), intercept=FALSE)

  expect_equivalent(coef(glm_fit),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(ols_cost, sgs_cost, tol=1e-3)
})
