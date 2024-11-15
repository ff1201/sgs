test_that("test screening returns same output for SGS with l2 standardisation and intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  groups = rep(1:20,each=5)
  path_length = 10
  sgs_screen = fit_sgs(X=X,y=y, groups=groups, type="linear",alpha=0.95, path_length = 10,vFDR=0.1, gFDR=0.1,standardise="l2",intercept=TRUE,screen=TRUE)
  sgs_no_screen = fit_sgs(X=X,y=y, groups=groups, type="linear",alpha=0.95, path_length = 10,vFDR=0.1, gFDR=0.1,standardise="l2",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(sgs_screen$beta),
    as.matrix(sgs_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for gSLOPE with l2 standardisation and intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  groups = rep(1:20,each=5)
  path_length = 10
  gslope_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10,gFDR=0.1,standardise="l2",intercept=TRUE,screen=TRUE)
  gslope_no_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10,gFDR=0.1,standardise="l2",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(gslope_screen$beta),
    as.matrix(gslope_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGS with l2 standardisation and no intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgs_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="l2",intercept=FALSE,screen=TRUE)
  sgs_no_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(sgs_screen$beta),
    as.matrix(sgs_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for gSLOPE with l2 standardisation and no intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  gslope_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10, gFDR=0.1,standardise="l2",intercept=FALSE,screen=TRUE)
  gslope_no_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10, gFDR=0.1,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(gslope_screen$beta),
    as.matrix(gslope_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGS with no standardisation and an intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgs_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="none",intercept=TRUE,screen=TRUE)
  sgs_no_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="none",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(sgs_screen$beta),
    as.matrix(sgs_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for gSLOPE with no standardisation and an intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  gslope_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10, gFDR=0.1,standardise="none",intercept=TRUE,screen=TRUE)
  gslope_no_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10, gFDR=0.1,standardise="none",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(gslope_screen$beta),
    as.matrix(gslope_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGS with no standardisation and no intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgs_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="none",intercept=FALSE,screen=TRUE)
  sgs_no_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="none",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(sgs_screen$beta),
    as.matrix(sgs_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for gSLOPE with no standardisation and no intercept", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  gslope_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10, gFDR=0.1,standardise="none",intercept=FALSE,screen=TRUE)
  gslope_no_screen = fit_gslope(X=X,y=y, groups=groups, type="linear", path_length = 10, gFDR=0.1,standardise="none",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(gslope_screen$beta),
    as.matrix(gslope_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGS with alpha = 0.05", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  path_length = 10
  groups = rep(1:20,each=5)
  sgs_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.05, vFDR=0.1, gFDR=0.1,standardise="l2",intercept=TRUE,screen=TRUE)
  sgs_no_screen = fit_sgs(X=X,y=y, groups=groups, type="linear", path_length = 10, alpha=0.05, vFDR=0.1, gFDR=0.1,standardise="l2",intercept=TRUE,screen=FALSE)

  expect_equivalent(as.matrix(sgs_screen$beta),
    as.matrix(sgs_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for SGS with logistic regression", {
  n = 50
  p = 100
  X = MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  path_length = 10
  groups = rep(1:20,each=5)
  sgs_screen = fit_sgs(X=X,y=y, groups=groups, type="logistic", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="l2",intercept=FALSE,screen=TRUE)
  sgs_no_screen = fit_sgs(X=X,y=y, groups=groups, type="logistic", path_length = 10, alpha=0.95, vFDR=0.1, gFDR=0.1,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(sgs_screen$beta),
    as.matrix(sgs_no_screen$beta),
    tol = 1e-3
  )
})

test_that("test screening returns same output for gSLOPE with logistic regression", {
  n = 50
  p = 100
  X = MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=diag(1,p))
  y = 1/(1+exp(-(X %*%rnorm(p,mean=0,sd=sqrt(10)) + rnorm(n,mean=0,sd=4))))
  y = ifelse(y>0.5,1,0)
  path_length = 10
  groups = rep(1:20,each=5)
  gslope_screen = fit_gslope(X=X,y=y, groups=groups, type="logistic", path_length = 10, gFDR=0.1,standardise="l2",intercept=FALSE,screen=TRUE)
  gslope_no_screen = fit_gslope(X=X,y=y, groups=groups, type="logistic", path_length = 10, gFDR=0.1,standardise="l2",intercept=FALSE,screen=FALSE)

  expect_equivalent(as.matrix(gslope_screen$beta),
    as.matrix(gslope_no_screen$beta),
    tol = 1e-3
  )
})