% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.sgs}
\alias{plot.sgs}
\title{Plot models of the following object types: \code{"sgs"}, \code{"sgs_cv"}, \code{"gslope"}, \code{"gslope_cv"}.}
\usage{
\method{plot}{sgs}(x, how_many = 10, ...)
}
\arguments{
\item{x}{Object of one of the following classes: \code{"sgs"}, \code{"sgs_cv"}, \code{"gslope"}, \code{"gslope_cv"}.}

\item{how_many}{Defines how many predictors to plot. Plots the predictors in decreasing order of largest absolute value.}

\item{...}{further arguments passed to base function.}
}
\value{
A list containing:
\item{response}{The predicted response. In the logistic case, this represents the predicted class probabilities.}
\item{class}{The predicted class assignments. Only returned if type = "logistic" in the model object.}
}
\description{
Plots the pathwise solution of a cross-validation fit, from a call to one of the following: \code{\link[=fit_sgs]{fit_sgs()}}, \code{\link[=fit_sgs_cv]{fit_sgs_cv()}}, \code{\link[=fit_gslope]{fit_gslope()}}, \code{\link[=fit_gslope_cv]{fit_gslope_cv()}}, \code{\link[=fit_sgo]{fit_sgo()}}, \code{\link[=fit_sgo_cv]{fit_sgo_cv()}}, \code{\link[=fit_goscar]{fit_goscar()}}, \code{\link[=fit_goscar_cv]{fit_goscar_cv()}}.
}
\examples{
# specify a grouping structure
groups = c(1,1,2,2,3)
# generate data
data =  gen_toy_data(p=5, n=4, groups = groups, seed_id=3,signal_mean=20,group_sparsity=1)
# run SGS 
model = fit_sgs(X = data$X, y = data$y, groups=groups, type = "linear", 
path_length = 20, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, 
min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=FALSE)
plot(model, how_many = 10)
}
\seealso{
\code{\link[=fit_sgs]{fit_sgs()}}, \code{\link[=fit_sgs_cv]{fit_sgs_cv()}}, \code{\link[=fit_gslope]{fit_gslope()}}, \code{\link[=fit_gslope_cv]{fit_gslope_cv()}}, \code{\link[=fit_sgo]{fit_sgo()}}, \code{\link[=fit_sgo_cv]{fit_sgo_cv()}}, \code{\link[=fit_goscar]{fit_goscar()}}, \code{\link[=fit_goscar_cv]{fit_goscar_cv()}}

Other SGS-methods: 
\code{\link{as_sgs}()},
\code{\link{coef.sgs}()},
\code{\link{fit_sgo}()},
\code{\link{fit_sgo_cv}()},
\code{\link{fit_sgs}()},
\code{\link{fit_sgs_cv}()},
\code{\link{predict.sgs}()},
\code{\link{print.sgs}()},
\code{\link{scaled_sgs}()}

Other gSLOPE-methods: 
\code{\link{coef.sgs}()},
\code{\link{fit_goscar}()},
\code{\link{fit_goscar_cv}()},
\code{\link{fit_gslope}()},
\code{\link{fit_gslope_cv}()},
\code{\link{predict.sgs}()},
\code{\link{print.sgs}()}
}
\concept{SGS-methods}
\concept{gSLOPE-methods}
