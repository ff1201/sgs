% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_goscar_cv.R
\name{fit_goscar_cv}
\alias{fit_goscar_cv}
\title{Fit a gOSCAR model using k-fold cross-validation.}
\usage{
fit_goscar_cv(
  X,
  y,
  groups,
  type = "linear",
  lambda = "path",
  path_length = 20,
  min_frac = 0.05,
  nfolds = 10,
  backtracking = 0.7,
  max_iter = 5000,
  max_iter_backtracking = 100,
  tol = 1e-05,
  standardise = "l2",
  intercept = TRUE,
  error_criteria = "mse",
  screen = TRUE,
  verbose = FALSE,
  w_weights = NULL,
  warm_start = NULL
)
}
\arguments{
\item{X}{Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).}

\item{y}{Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.}

\item{groups}{A grouping structure for the input data. Should take the form of a vector of group indices.}

\item{type}{The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.}

\item{lambda}{The regularisation parameter. Defines the level of sparsity in the model. A higher value leads to sparser models:
\itemize{
\item \code{"path"} computes a path of regularisation parameters of length \code{"path_length"}. The path will begin just above the value at which the first predictor enters the model and will terminate at the value determined by \code{"min_frac"}.
\item User-specified single value or sequence. Internal scaling is applied based on the type of standardisation. The returned \code{"lambda"} value will be the original unscaled value(s).
}}

\item{path_length}{The number of \eqn{\lambda} values to fit the model for. If \code{"lambda"} is user-specified, this is ignored.}

\item{min_frac}{Smallest value of \eqn{\lambda} as a fraction of the maximum value. That is, the final \eqn{\lambda} will be \code{"min_frac"} of the first \eqn{\lambda} value.}

\item{nfolds}{The number of folds to use in cross-validation.}

\item{backtracking}{The backtracking parameter, \eqn{\tau}, as defined in Pedregosa and Gidel (2018).}

\item{max_iter}{Maximum number of ATOS iterations to perform.}

\item{max_iter_backtracking}{Maximum number of backtracking line search iterations to perform per global iteration.}

\item{tol}{Convergence tolerance for the stopping criteria.}

\item{standardise}{Type of standardisation to perform on \code{X}:
\itemize{
\item \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one.
\item \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one.
\item \code{"sd"} standardises the input data to have standard deviation of one.
\item \code{"none"} no standardisation applied.
}}

\item{intercept}{Logical flag for whether to fit an intercept.}

\item{error_criteria}{The criteria used to discriminate between models along the path. Supported values are: \code{"mse"} (mean squared error) and \code{"mae"} (mean absolute error).}

\item{screen}{Logical flag for whether to apply screening rules (see Feser and Evangelou (2024)). Screening discards irrelevant groups before fitting, greatly improving speed.}

\item{verbose}{Logical flag for whether to print fitting information.}

\item{w_weights}{Optional vector for the group penalty weights. Overrides the OSCAR penalties when specified. When entering custom weights, these are multiplied internally by \eqn{\lambda}. To void this behaviour, set \eqn{\lambda = 1}.}

\item{warm_start}{Optional list for implementing warm starts. These values are used as initial values in the fitting algorithm. Need to supply \code{"x"} and \code{"u"} in the form \code{"list(warm_x, warm_u)"}. Not recommended for use with a path or CV fit as start from the null model by design.}
}
\value{
A list containing:
\item{errors}{A table containing fitting information about the models on the path.}
\item{all_models}{Fitting information for all models fit on the path, which is a \code{"gslope"} object type.}
\item{fit}{The 1se chosen model, which is a \code{"gslope"} object type.}
\item{best_lambda}{The value of \eqn{\lambda} which generated the chosen model.}
\item{best_lambda_id}{The path index for the chosen model.}
}
\description{
Function to fit a pathwise solution of group OSCAR (gOSCAR) models using k-fold cross-validation. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
}
\details{
Fits gOSCAR models under a pathwise solution using adaptive three operator splitting (ATOS), picking the 1se model as optimum. Warm starts are implemented.
}
\examples{
# specify a grouping structure
groups = c(1,1,1,2,2,3,3,3,4,4)
# generate data
data =  gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
# run gOSCAR with cross-validation
cv_model = fit_goscar_cv(X = data$X, y = data$y, groups=groups, type = "linear", path_length = 5, 
nfolds=5, min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=TRUE)
}
\references{
Bao, R., Gu B., Huang, H. (2020). \emph{Fast OSCAR and OWL Regression via Safe Screening Rules}, \url{https://proceedings.mlr.press/v119/bao20b}

Feser, F., Evangelou, M. (2024). \emph{Strong screening rules for group-based SLOPE models}, \url{https://arxiv.org/abs/2405.15357}
}
\seealso{
\code{\link[=fit_goscar]{fit_goscar()}}

Other gSLOPE-methods: 
\code{\link{coef.sgs}()},
\code{\link{fit_goscar}()},
\code{\link{fit_gslope}()},
\code{\link{fit_gslope_cv}()},
\code{\link{plot.sgs}()},
\code{\link{predict.sgs}()},
\code{\link{print.sgs}()}

Other model-selection: 
\code{\link{as_sgs}()},
\code{\link{fit_gslope_cv}()},
\code{\link{fit_sgo_cv}()},
\code{\link{fit_sgs_cv}()},
\code{\link{scaled_sgs}()}
}
\concept{gSLOPE-methods}
\concept{model-selection}
