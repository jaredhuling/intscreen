
#' filters top interactions using correlation matrix
#' @description This function filters top interactions using correlation matrix
#' @param cormat correlation matrix returned by \code{\link[intscreen]{compute_cors}}
#' @param nints integer number of top interactions to screen
#' @seealso \code{\link[intscreen]{compute_cors}}, \code{\link[intscreen]{construct_ints}}, and \code{\link[intscreen]{intscreen}}
#' @export
filter_top <- function(cormat, nints = 50)
{
    corvals <- cormat[upper.tri(cormat, diag = TRUE)]
    corvals <- corvals[(!is.na(corvals))]
    corvals <- abs(corvals[(!is.nan(corvals))])

    num_ints <- length(corvals)
    if (any(nints > num_ints))
    {
        stop("nints is larger than the number of total interactions. choose a smaller number")
    }

    quantval <- quantile(corvals, probs = 1 - nints / num_ints)

    top_ints <- which(abs(cormat) >= quantval - 1e-12, arr.ind = TRUE)

    top_ints
}


#' computes unscaled correlation matrix
#' @description This function computes an unscaled correlation matrix
#' @param x matrix of predictors
#' @param y vector of observations
#' @param modifier effect modifier
#' @param subset1 first subset of variable indices
#' @param subset2 second subset of variable indices
#' @seealso \code{\link[intscreen]{filter_top}}, \code{\link[intscreen]{construct_ints}}, and \code{\link[intscreen]{intscreen}}
#' @export
compute_cors <- function(x, y, modifier = NULL, subset1 = NULL, subset2 = NULL)
{
    if (is.null(subset1) & is.null(subset2))
    {
        if (is.null(modifier))
        {
            return(compute_cors_cpp(X = x, Y = scale(y)))
        } else
        {
            return(compute_cors_mod_cpp(X = x, Y = scale(y), mod = modifier))
        }
    } else if (is.null(subset1) & !is.null(subset2))
    {
        subset1 <- as.integer(subset2)
        subset2 <- as.integer(1:ncol(x))
        if (is.null(modifier))
        {
            return(compute_cors_subset_cpp(X = x, Y = scale(y),
                                           idx1 = subset1, idx2 = subset2))
        } else
        {
            return(compute_cors_subset_mod_cpp(X = x, Y = scale(y), mod = modifier,
                                               idx1 = subset1, idx2 = subset2))
        }
    } else if (!is.null(subset1) & is.null(subset2))
    {
        subset2 <- as.integer(1:ncol(x))
        subset1 <- as.integer(subset1)
        if (is.null(modifier))
        {
            return(compute_cors_subset_cpp(X = x, Y = scale(y),
                                           idx1 = subset1, idx2 = subset2))
        } else
        {
            return(compute_cors_subset_mod_cpp(X = x, Y = scale(y), mod = modifier,
                                               idx1 = subset1, idx2 = subset2))
        }
    } else (!is.null(subset1) & !is.null(subset2))
    {
        subset1 <- as.integer(subset1)
        subset2 <- as.integer(subset2)
        if (is.null(modifier))
        {
            return(compute_cors_subset_cpp(X = x, Y = scale(y),
                                           idx1 = subset1, idx2 = subset2))
        } else
        {
            return(compute_cors_subset_mod_cpp(X = x, Y = scale(y), mod = modifier,
                                               idx1 = subset1, idx2 = subset2))
        }
    }
}



#' constructs selected interactions
#' @description This function constructs interaction terms
#' @param x matrix of predictors
#' @param top_ints \code{nints} x 2 integer matrix of interaction terms
#' @param modifier effect modifier
#' @seealso \code{\link[intscreen]{filter_top}}, \code{\link[intscreen]{compute_cors}}, and \code{\link[intscreen]{intscreen}}
#' @export
#' @examples
#'
#' x <- matrix(rnorm(1000 * 500), ncol = 500)
#' y <- rnorm(1000)
#'
#' cormat <- compute_cors(x, y)
#'
#' int_idx <- filter_top(cormat, nints = 50)
#'
#' int_mat <- construct_ints(x, ints)
construct_ints <- function(x, top_ints, modifier = NULL)
{
    cnames <- colnames(x)

    stopifnot(is.matrix(top_ints))

    stopifnot(ncol(top_ints) == 2)

    if (mode(top_ints) != "integer")
    {
        mode(top_ints) <- "integer"
    }

    if (is.null(cnames))
    {
        cnames <- paste0("X", 1:ncol(x))
    }

    intnames <- paste(cnames[top_ints[,1]], cnames[top_ints[,2]], sep = ":")

    if (is.null(modifier))
    {
        int_mat <- construct_ints_cpp(X = x, whichints = top_ints)
    } else
    {
        int_mat <- modifier * construct_ints_cpp(X = x, whichints = top_ints)
    }

    colnames(int_mat) <- intnames

    int_mat
}



#' all computations for interaction screening
#' @description This function does all computations for interaction screening
#' @param x matrix of predictors
#' @param y vector of observations
#' @param nints integer number of top interactions to screen
#' @param heredity either \code{"weak"}, \code{"strong"}, or \code{"none"}
#' @param modifier effect modifier
#' @export
#' @examples
#'
#' x <- matrix(rnorm(1000 * 500), ncol = 500)
#' y <- rnorm(1000)
#'
#' ints <- intscreen(x, y, nints = 12)
#'
#' str(ints)
intscreen <- function(x, y, nints = 10,
                      heredity = c("none", "weak", "strong"),
                      modifier = NULL)
{
    heredity <- match.arg(heredity)

    if (heredity == "strong")
    {
        corxy     <- drop(cor(y, x))
        nints_1 <- min(nints, ncol(x))
        top_nints <- order(abs(corxy), decreasing = TRUE)[1:nints_1]

        cormat    <- compute_cors(x[,top_nints,drop=FALSE], y, modifier)

        int_idx_top <- filter_top(cormat, nints = nints)

        int_idx <- cbind(top_nints[int_idx_top[,1]], top_nints[int_idx_top[,2]])

        int_mat <- construct_ints(x, int_idx, modifier = modifier)

    } else if (heredity == "weak")
    {
        corxy     <- drop(cor(y, x))

        nints_1 <- min(nints, ncol(x))
        top_nints <- order(abs(corxy), decreasing = TRUE)[1:nints_1]

        non_top_nints  <- (1:ncol(x))[-top_nints]
        reordered_cols <- c(top_nints, non_top_nints)

        cormat    <- compute_cors(x, y, modifier,
                                  subset1 = top_nints,
                                  subset2 = reordered_cols)

        int_idx_top <- filter_top(cormat, nints = nints)

        int_idx <- cbind(top_nints[int_idx_top[,1]],
                         reordered_cols[int_idx_top[,2]])

        int_mat <- construct_ints(x, int_idx, modifier = modifier)
    } else if (heredity == "none")
    {
        cormat <- compute_cors(x, y, modifier)

        int_idx <- filter_top(cormat, nints = nints)

        int_mat <- construct_ints(x, int_idx, modifier = modifier)
    }

    list(int_mat = int_mat, int_idx = int_idx)
}


#' Cross validation interaction screening
#' @description This function implements CV interaction screening
#' @param x matrix of predictors
#' @param y vector of observations
#' @param nints integer number of top interactions to screen
#' @param nsplits integer number of cross validation splits to run. defaults to 10
#' @param train.frac fraction of data used for each split. defaults to 0.75. Only
#' used for \code{resample.type = "cv"}
#' @param fraction.in.thresh fraction of times across the \code{nsplits} CV splits each
#' interaction is required in the top \code{k} interactions in order to be selected
#' @param resample.type either \code{"cv"} for cross validation or \code{"bootstrap"} for
#' bootstrap approach
#' @param heredity either \code{"weak"}, \code{"strong"}, or \code{"none"}
#' @param verbose logical value, whether to print progress of the CV splitting
#' @param modifier effect modifier
#' @export
#' @examples
#'
#' set.seed(1)
#' x <- matrix(rnorm(150 * 500), ncol = 500)
#' y <- rnorm(150) + x[,1] * x[,2] * 0.5 - x[,3] * x[,4] * 0.5
#'
#' ## require that each interaction be in the top 200 ints 5% of the 10 splits
#' ints <- cv_intscreen(x, y, nints = 200, nsplits = 10, fraction.in.thresh = 0.75)
#'
#' ints$int_idx
#'
#' ## require that each interaction be in the top 50 ints 50% of the 10 splits
#' ints <- cv_intscreen(x, y, nints = 50, nsplits = 10, fraction.in.thresh = 0.5)
#'
#' ints$int_idx
cv_intscreen <- function(x, y, nints = 100, nsplits = 10, train.frac = 0.75, fraction.in.thresh = 1,
                         resample.type = c("bootstrap", "cv"),
                         heredity = c("none", "weak", "strong"),
                         verbose = FALSE, modifier = NULL)
{
    heredity     <- match.arg(heredity)
    resampletype <- match.arg(resample.type)

    intlist <- melist <- vector(mode = "list", length = nsplits)

    n <- nrow(x)

    for (s in 1:nsplits)
    {
        if (resampletype == "cv")
        {
            s_idx     <- sample.int(n, floor(n * train.frac))
        } else
        {
            s_idx     <- sample.int(n, n, replace = TRUE)
        }

        corxy     <- drop(cor(y[s_idx], x[s_idx,,drop=FALSE]))
        nints_1   <- min(nints, ncol(x))
        top_nints <- order(abs(corxy), decreasing = TRUE)[1:nints_1]

        melist[[s]] <- top_nints

        if (heredity == "strong")
        {

            if (is.null(modifier))
            {
                cormat <- compute_cors(x[s_idx,top_nints,drop=FALSE], y[s_idx])
            } else
            {
                cormat <- compute_cors(x[s_idx,top_nints,drop=FALSE], y[s_idx], modifier[s_idx])
            }

            int_idx_top  <- filter_top(cormat, nints = nints)

            intlist[[s]] <- cbind(top_nints[int_idx_top[,1]], top_nints[int_idx_top[,2]])

        } else if (heredity == "weak")
        {

            non_top_nints  <- (1:ncol(x))[-top_nints]
            reordered_cols <- c(top_nints, non_top_nints)

            if (is.null(modifier))
            {
                cormat    <- compute_cors(x[s_idx,,drop=FALSE],
                                          y[s_idx],
                                          subset1 = top_nints,
                                          subset2 = reordered_cols)
            } else
            {
                cormat    <- compute_cors(x[s_idx,,drop=FALSE],
                                          y[s_idx],
                                          modifier[s_idx],
                                          subset1 = top_nints,
                                          subset2 = reordered_cols)
            }

            int_idx_top <- filter_top(cormat, nints = nints)

            intlist[[s]] <- cbind(top_nints[int_idx_top[,1]],
                                  reordered_cols[int_idx_top[,2]])
        } else if (heredity == "none")
        {
            if (is.null(modifier))
            {
                cormat  <- compute_cors(x[s_idx,,drop=FALSE], y[s_idx])
            } else
            {
                cormat  <- compute_cors(x[s_idx,,drop=FALSE], y[s_idx], modifier[s_idx])
            }

            intlist[[s]] <- filter_top(cormat, nints = nints)
        }


        if (verbose)
        {
            cat("Split number:", s, "\n")
        }
    }

    ## the first index (ie first column of matrix returned by filter_top)
    ## is always less than the second index, so we don't need to worry
    ## about ordering of interactions
    int_tall <- do.call(rbind, intlist)

    all_ints <- paste(int_tall[,1], int_tall[,2], sep = ":")

    tab <- table(all_ints)

    ints_all_in <- tab[tab >= fraction.in.thresh * nsplits]


    me_all <- unlist(melist)
    tab    <- table(me_all)

    me_all_in <- tab[tab >= fraction.in.thresh * nsplits]

    me_all_in <- as.integer(names(me_all_in))

    if (length(ints_all_in) == 0)
    {
        warning("No interactions found.")
        int_idx_final <- NULL
        int_mat <- NULL
    } else
    {
        ints_in_char <- names(ints_all_in)

        int_idx_final <- do.call(rbind, lapply(strsplit(ints_in_char, ":"), as.integer))

        colnames(int_idx_final) <- colnames(int_tall)
        int_mat <- construct_ints(x, int_idx_final, modifier = modifier)
    }

    if (length(me_all_in) == 0)
    {
        me_all_in <- NULL
    }

    list(int_mat = int_mat, int_idx = int_idx_final, main_effects = me_all_in)
}


