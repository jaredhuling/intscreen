

#' Cross validation interaction screening
#' @description This function implements CV interaction screening
#' @param x matrix of predictors
#' @param y vector of observations
#' @param k integer number of top interactions to screen
#' @param nsplits integer number of cross validation splits to run. defaults to 10
#' @param train.frac fraction of data used for each split. defaults to 0.75
#' @param fraction.in.thresh fraction of times across the \code{nsplits} CV splits each
#' interaction is required in the top \code{k} interactions in order to be selected
#' @param ... arguments to be passed to \code{\link[glmnet]{glmnet}}
#' @export
#' @examples
#'
#' set.seed(1)
#' x <- matrix(rnorm(100 * 350), ncol = 350)
#' y <- rnorm(100) + x[,1] * x[,2] * 0.5 - x[,3] * x[,4] * 0.5
#'
#'
#' ## require that each interaction be in the top 50 ints 100% of the 15 splits
#' intmod <- intnet(x, y, k = 50, nsplits = 15, fraction.in.thresh = 1)
#'
#' plot(intmod)
intnet <- function(x, y, k = 100, nsplits = 10, train.frac = 0.75, fraction.in.thresh = 1,
                   verbose = FALSE, ...)
{
    cnames <- colnames(x)

    if (is.null(cnames))
    {
        cnames <- paste0("X", 1:ncol(x))
        colnames(x) <- cnames
    }

    ints <- cv_intscreen(x = x, y = y, k = k, nsplits = nsplits, train.frac = train.frac,
                         fraction.in.thresh = fraction.in.thresh, verbose = FALSE)

    if (!is.null(ints$int_mat))
    {
        x_mm <- cbind(x, ints$int_mat)
        obj  <- glmnet(x = x_mm, y = y, ...)
    } else
    {
        obj  <- glmnet(x = x, y = y, ...)
    }
    obj$ints <- ints
    class(obj)[class(obj) == "glmnet"] <- "intnet"
    class(obj) <- c(class(obj), "glmnet")
    obj
}



#' Cross validation interaction screening
#' @description This function implements CV interaction screening
#' @param ... other arguments to be passed to \code{\link[intscreen]{intnet}}
#' @inheritParams glmnet::cv.glmnet
#' @export
#' @import glmnet
#' @examples
#'
#' set.seed(1)
#' x <- matrix(rnorm(100 * 350), ncol = 350)
#' y <- rnorm(100) + x[,1] * x[,2] * 0.5 - x[,3] * x[,4] * 0.5
#'
#'
#' ## require that each interaction be in the top 50 ints 100% of the 15 splits
#' intmod <- cv.intnet(x, y, k = 50, nsplits = 15, fraction.in.thresh = 1)
#'
#' plot(intmod)
#'
#' cfs <- as.matrix(predict(intmod, type = "coef", s = "lambda.min"))
#' cfs[cfs != 0,,drop=FALSE]
cv.intnet <- function (x, y, weights, offset = NULL, lambda = NULL, type.measure = c("mse",
                                                                                     "deviance", "class", "auc", "mae"), nfolds = 10, foldid,
                       grouped = TRUE, keep = FALSE, parallel = FALSE, ...)
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

    glmnet.object = intnet(x, y, weights = weights, offset = offset,
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
            intnet(x[!which, , drop = FALSE], y_sub, lambda = lambda,
                   offset = offset_sub, weights = weights[!which],
                   ...)
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
            outlist[[i]] = intnet(x[!which, , drop = FALSE],
                                  y_sub, lambda = lambda, offset = offset_sub,
                                  weights = weights[!which], ...)
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

        newx <- cbind(newx, intmat)
    }


    nfit=as.matrix(cbind2(1,newx)%*%nbeta)
    if(object$offset){
        if(missing(newoffset))stop("No newoffset provided for prediction, yet offset used in fit of glmnet",call.=FALSE)
        if(is.matrix(newoffset)&&inherits(object,"lognet")&&dim(newoffset)[[2]]==2)newoffset=newoffset[,2]
        nfit=nfit+array(newoffset,dim=dim(nfit))
    }
    nfit
}


