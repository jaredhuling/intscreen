% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intnet.R
\name{cv.intnet}
\alias{cv.intnet}
\title{Cross validation interaction screening}
\usage{
cv.intnet(x, y, weights, offset = NULL, lambda = NULL,
  type.measure = c("mse", "deviance", "class", "auc", "mae"),
  nfolds = 10, foldid, grouped = TRUE, keep = FALSE,
  parallel = FALSE, modifier = NULL, ...)
}
\arguments{
\item{x}{\code{x} matrix as in \code{glmnet}.}

\item{y}{response \code{y} as in \code{glmnet}.}

\item{weights}{Observation weights; defaults to 1 per observation}

\item{offset}{Offset vector (matrix) as in \code{glmnet}}

\item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{glmnet} chooses its own sequence}

\item{type.measure}{loss to use for cross-validation. Currently five
    options, not all available for all models. The default is \code{type.measure="deviance"}, which uses
    squared-error for gaussian models (a.k.a \code{type.measure="mse"} there), deviance
    for logistic and poisson
    regression, and partial-likelihood for the Cox
    model. \code{type.measure="class"} applies to binomial and multinomial logistic regression only,
    and gives misclassification error.  \code{type.measure="auc"} is for
    two-class logistic regression only, and gives area under the ROC
    curve. \code{type.measure="mse"} or \code{type.measure="mae"} (mean absolute error)
    can be used by all models except the \code{"cox"}; they measure the
    deviation from the fitted mean to the response.}

\item{nfolds}{number of folds - default is 10. Although \code{nfolds}
    can be as large as the sample size (leave-one-out CV), it is not
    recommended for large datasets. Smallest value allowable is \code{nfolds=3}}

\item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying what fold each observation is in. If supplied,
    \code{nfold} can be missing.}

\item{grouped}{This is an experimental argument, with default
      \code{TRUE}, and can be ignored by most users. For all models
      except the \code{"cox"}, this refers to computing \code{nfolds}
      separate statistics, and then using their mean and estimated
      standard error to describe the CV curve. If \code{grouped=FALSE},
      an error matrix is built up at the observation level from the predictions
      from the \code{nfold} fits, and then summarized (does not apply to
      \code{type.measure="auc"}). For the \code{"cox"} family,
      \code{grouped=TRUE} obtains the CV partial likelihood for the Kth
      fold by \emph{subtraction}; by subtracting the log partial
      likelihood evaluated on the full dataset from that evaluated on
      the  on the (K-1)/K dataset. This makes more efficient use of risk
      sets. With \code{grouped=FALSE} the log partial likelihood is
      computed only on the Kth fold}

\item{keep}{If \code{keep=TRUE}, a \emph{prevalidated} array is
      returned containing fitted values for each observation and each
      value of \code{lambda}. This means these fits are computed with
      this observation and the rest of its fold omitted. The
      \code{folid} vector is also returned. Default is {keep=FALSE}}

\item{parallel}{If \code{TRUE}, use parallel \code{foreach} to fit each fold.
      Must register parallel before hand, such as \code{doMC} or others.
    See the example below.}

\item{modifier}{effect modifier}

\item{...}{other arguments to be passed to \code{\link[intscreen]{intnet}}}
}
\description{
This function implements CV interaction screening
}
\examples{

library(intscreen)

set.seed(1)
x <- matrix(rnorm(100 * 350), ncol = 350)
y <- rnorm(100) + x[,1] * x[,2] * 0.5 - x[,3] * x[,4] * 0.5


## require that each interaction be in the top 50 ints 100\% of the 15 splits
intmod <- cv.intnet(x, y, nints = 50, nsplits = 15, fraction.in.thresh = 1)

plot(intmod)

cfs <- as.matrix(predict(intmod, type = "coef", s = "lambda.min"))
cfs[cfs != 0,,drop=FALSE]
}
