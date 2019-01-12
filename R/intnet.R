

#' Cross validation interaction screening
#' @description This function implements CV interaction screening
#' @param x matrix of predictors
#' @param y vector of observations
#' @param which.cols integer vector indicating which columns of \code{x} to check for interactions
#' @param nints integer number of top interactions to screen
#' @param resampletype either \code{"cv"} for cross validation or \code{"bootstrap"} for
#' bootstrap approach
#' @param heredity either \code{"weak"}, \code{"strong"}, or \code{"none"}
#' @param nsplits integer number of cross validation splits to run. defaults to 10
#' @param train.frac fraction of data used for each split. defaults to 0.75
#' @param fraction.in.thresh fraction of times across the \code{nsplits} CV splits each
#' interaction is required in the top \code{nints} interactions in order to be selected
#' @param modifier effect modifier
#' @param ... arguments to be passed to \code{\link[glmnet]{glmnet}}
#' @export
#' @examples
#' library(intnet)
#'
#' set.seed(1)
#' x <- matrix(rnorm(100 * 350), ncol = 350)
#' y <- rnorm(100) + x[,1] * x[,2] * 0.5 - x[,3] * x[,4] * 0.5
#'
#'
#' ## require that each interaction be in the top 50 ints 100% of the 15 splits
#' intmod <- intnet(x, y, nints = 50, nsplits = 15, fraction.in.thresh = 1)
#'
#' plot(intmod)
intnet <- function(x, y, which.cols = 1:ncol(x),
                   nints = 100,
                   heredity = c("none", "weak", "strong"),
                   resample.type = c("bootstrap", "cv"),
                   nsplits = 10,
                   train.frac = 0.75, fraction.in.thresh = 1,
                   verbose = FALSE,
                   modifier = NULL,
                   ...)
{
    heredity     <- match.arg(heredity)
    resampletype <- match.arg(resample.type)

    cnames   <- colnames(x)

    if (is.null(cnames))
    {
        cnames <- paste0("X", 1:ncol(x))
        colnames(x) <- cnames
    }

    if (length(which.cols) < 2)
    {
        stop("which.cols must be a vector of at least length 2")
    } else
    {
        which.cols <- as.integer(which.cols)
    }

    if (max(which.cols) > ncol(x))
    {
        stop("which.cols should contain values no larger than the number of columns in x")
    }

    if (min(which.cols) < 1)
    {
        stop("which.cols should contain values no smaller than 1")
    }

    if (is.null(modifier))
    {
        modifier <- rep(1, nrow(x))
    }

    ints <- cv_intscreen(x = x[,which.cols,drop=FALSE], y = y, nints = nints,
                         nsplits = nsplits, train.frac = train.frac,
                         fraction.in.thresh = fraction.in.thresh,
                         resampletype = resampletype,
                         heredity = heredity,
                         verbose = FALSE,
                         modifier = modifier)

    if (is.null(modifier))
    {
        if (!is.null(ints$main_effects))
        {
            x <- x[,ints$main_effects,drop=FALSE]
        }

        if (!is.null(ints$int_mat))
        {
            x_mm <- cbind(x, ints$int_mat)
            obj  <- glmnet(x = x_mm, y = y, ...)
        } else
        {
            obj  <- glmnet(x = x, y = y, ...)
        }
    } else
    {

        if (!is.null(ints$main_effects))
        {
            x <- x[,ints$main_effects,drop=FALSE]
        }

        if (!is.null(ints$int_mat))
        {
            x_mm <- modifier * cbind(x, ints$int_mat)
            obj  <- glmnet(x = x_mm, y = y, ...)
        } else
        {
            obj  <- glmnet(x = modifier * x, y = y, ...)
        }
    }

    obj$ints <- ints
    class(obj)[class(obj) == "glmnet"] <- "intnet"
    class(obj) <- c(class(obj), "glmnet")
    obj
}

## the following code is modified from glmnet:

#' Cross validation interaction screening
#' @description This function implements CV interaction screening
#' @param ... other arguments to be passed to \code{\link[intscreen]{intnet}}
#' @param modifier effect modifier
#' @inheritParams glmnet::cv.glmnet
#' @export
#' @import glmnet
#' @examples
#'
#' library(intscreen)
#'
#' set.seed(1)
#' x <- matrix(rnorm(100 * 350), ncol = 350)
#' y <- rnorm(100) + x[,1] * x[,2] * 0.5 - x[,3] * x[,4] * 0.5
#'
#'
#' ## require that each interaction be in the top 50 ints 100% of the 15 splits
#' intmod <- cv.intnet(x, y, nints = 50, nsplits = 15, fraction.in.thresh = 1)
#'
#' plot(intmod)
#'
#' cfs <- as.matrix(predict(intmod, type = "coef", s = "lambda.min"))
#' cfs[cfs != 0,,drop=FALSE]
cv.intnet <- function (x, y, weights, offset = NULL, lambda = NULL,
                       type.measure = c("mse", "deviance", "class", "auc", "mae"),
                       nfolds = 10, foldid,
                       grouped = TRUE, keep = FALSE, parallel = FALSE, modifier = NULL, ...)
{
    if (missing(type.measure))
        type.measure = "default"
    else type.measure = match.arg(type.measure)
    if (!is.null(lambda) && length(lambda) < 2)
        stop("Need more than one value of lambda for cv.glmnet")
    N = nrow(x)
    if (missing(weights))
        weights = rep(1, N)
    else weights = as.double(weights)
    y = drop(y)
    glmnet.call = match.call(expand.dots = TRUE)
    which = match(c("type.measure", "nfolds", "foldid", "grouped",
                    "keep"), names(glmnet.call), FALSE)
    if (any(which))
        glmnet.call = glmnet.call[-which]
    glmnet.call[[1]] = as.name("glmnet")


    dots <- list(...)

    intscr_args <- c("k", "nsplits", "train.frac", "fraction.in.thresh", "verbose")

    if (any(intscr_args %in% names(dots)))
    {
        instr_dots <- dots[intscr_args %in% names(dots)]
    }

    glmnet.object = intnet(x, y, weights = weights,
                           offset = offset,
                           modifier = modifier,
                           lambda = lambda, ...)
    glmnet.object$call = glmnet.call
    subclass=class(glmnet.object)[[1]]
    type.measure=cvtype(type.measure,subclass)
    is.offset = glmnet.object$offset
    ###Next line is commented out so each call generates its own lambda sequence
    # lambda=glmnet.object$lambda
    if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
        nz = predict(glmnet.object, type = "nonzero")
        nz = sapply(nz, function(x) sapply(x, length))
        nz = ceiling(apply(nz, 1, median))
    }
    else nz = sapply(predict(glmnet.object, type = "nonzero"),
                     length)
    if (missing(foldid))
        foldid = sample(rep(seq(nfolds), length = N))
    else nfolds = max(foldid)
    if (nfolds < 3)
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist = as.list(seq(nfolds))
    if (parallel) {
        #  if (parallel && require(foreach)) {
        outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
        {
            which = foldid == i
            #      if (is.matrix(y))
            if (length(dim(y))>1)
                y_sub = y[!which, ]
            else y_sub = y[!which]
            if (is.offset)
                offset_sub = as.matrix(offset)[!which, ]
            else offset_sub = NULL
            if (is.null(modifier))
            {
                return(intnet(x[!which, , drop = FALSE], y_sub, lambda = lambda,
                              offset = offset_sub, weights = weights[!which],
                              ...))
            } else
            {
                return(intnet(x[!which, , drop = FALSE], y_sub, lambda = lambda,
                              modifier = modifier[!which],
                              offset = offset_sub, weights = weights[!which],
                              ...))
            }
        }
    }
    else {
        for (i in seq(nfolds)) {
            which = foldid == i
            if (is.matrix(y))
                y_sub = y[!which, ]
            else y_sub = y[!which]
            if (is.offset)
                offset_sub = as.matrix(offset)[!which, ]
            else offset_sub = NULL
            if (is.null(modifier))
            {
                outlist[[i]] = intnet(x[!which, , drop = FALSE],
                                      y_sub, lambda = lambda, offset = offset_sub,
                                      weights = weights[!which], ...)
            } else
            {
                outlist[[i]] = intnet(x[!which, , drop = FALSE],
                                      y_sub, lambda = lambda, offset = offset_sub,
                                      modifier = modifier[!which],
                                      weights = weights[!which], ...)
            }
        }
    }
    fun = paste("cv", subclass, sep = ".")
    lambda = glmnet.object$lambda
    cvstuff = do.call(fun, list(outlist, lambda, x, y, weights,
                                offset, foldid, type.measure, grouped, keep))
    cvm = cvstuff$cvm
    cvsd = cvstuff$cvsd
    nas=is.na(cvsd)
    if(any(nas)){
        lambda=lambda[!nas]
        cvm=cvm[!nas]
        cvsd=cvsd[!nas]
        nz=nz[!nas]
    }
    cvname = names(cvstuff$type.measure)
    names(cvname)=cvstuff$type.measure# to be compatible with earlier version; silly, I know
    out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
                   cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, glmnet.fit = glmnet.object)
    if (keep)
        out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
    lamin=if(cvname=="AUC")getmin(lambda,-cvm,cvsd)
    else getmin(lambda, cvm, cvsd)
    obj = c(out, as.list(lamin))
    class(obj) = "cv.glmnet"
    obj
}

#' prediction for intnet objects
#' @inheritParams glmnet::predict.glmnet
#' @importFrom stats approx predict
#' @method predict intnet
#' @export
predict.intnet=function(object,newx,s=NULL,type=c("link","response","coefficients","nonzero","class"),exact=FALSE,newoffset,...){
    type=match.arg(type)
    if(missing(newx)){
        if(!match(type,c("coefficients","nonzero"),FALSE))stop("You need to supply a value for 'newx'")
    }
    if(exact&&(!is.null(s))){
        ###we augment the lambda sequence with the new values, if they are different,and refit the model using update
        lambda=object$lambda
        which=match(s,lambda,FALSE)
        if(!all(which>0)){
            lambda=unique(rev(sort(c(s,lambda))))
            check_dots(object,...)# This fails if you don't supply the crucial arguments
            object=update(object,lambda=lambda,...)
        }
    }
    a0=t(as.matrix(object$a0))
    rownames(a0)="(Intercept)"
    nbeta=methods::rbind2(a0,object$beta)#was rbind2
    if(!is.null(s)){
        vnames=dimnames(nbeta)[[1]]
        dimnames(nbeta)=list(NULL,NULL)
        lambda=object$lambda
        lamlist=lambda.interp(lambda,s)

        nbeta=nbeta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac) +nbeta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
        dimnames(nbeta)=list(vnames,paste(seq(along=s)))
    }
    if(type=="coefficients")return(nbeta)
    if(type=="nonzero")return(nonzeroCoef(nbeta[-1,,drop=FALSE],bystep=TRUE))
    ###Check on newx
    if(inherits(newx, "sparseMatrix"))newx=as(newx,"dgeMatrix")

    if (!is.null(object$ints$int_idx))
    {
        intmat <- construct_ints(newx, object$ints$int_idx)

        if (!is.null(object$ints$main_effects))
        {
            newx <- cbind(newx[,object$ints$main_effects,drop=FALSE], intmat)
        } else
        {
            newx <- cbind(newx, intmat)
        }
    } else
    {
        if (!is.null(object$ints$main_effects))
        {
            newx <- newx[,object$ints$main_effects,drop=FALSE]
        }
    }

    nfit=as.matrix(cbind2(1,newx)%*%nbeta)
    if(object$offset){
        if(missing(newoffset))stop("No newoffset provided for prediction, yet offset used in fit of glmnet",call.=FALSE)
        if(is.matrix(newoffset)&&inherits(object,"lognet")&&dim(newoffset)[[2]]==2)newoffset=newoffset[,2]
        nfit=nfit+array(newoffset,dim=dim(nfit))
    }
    nfit
}



