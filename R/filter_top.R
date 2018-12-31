
#' filters top interactions using correlation matrix
#' @description This function filters top interactions using correlation matrix
#' @param cormat correlation matrix returned by \code{\link[intscreen]{compute_cors}}
#' @param k integer number of top interactions to screen
#' @seealso \code{\link[intscreen]{compute_cors}}, \code{\link[intscreen]{construct_ints}}, and \code{\link[intscreen]{intscreen}}
#' @export
filter_top <- function(cormat, k = 50)
{
    corvals <- cormat[upper.tri(cormat)]
    corvals <- corvals[(!is.na(corvals))]
    corvals <- abs(corvals[(!is.nan(corvals))])

    num_ints <- length(corvals)
    if (any(k > num_ints))
    {
        stop("k is larger than the number of total interactions. choose a smaller number")
    }

    quantval <- quantile(corvals, probs = 1 - k / num_ints)

    top_ints <- which(abs(cormat) >= quantval - 1e-12, arr.ind = TRUE)

    top_ints
}



#' computes unscaled correlation matrix
#' @description This function computes an unscaled correlation matrix
#' @param x matrix of predictors
#' @param y vector of observations
#' @seealso \code{\link[intscreen]{filter_top}}, \code{\link[intscreen]{construct_ints}}, and \code{\link[intscreen]{intscreen}}
#' @export
compute_cors <- function(x, y)
{
    compute_cors_cpp(X = x, Y = scale(y))
}


#' constructs selected interactions
#' @description This function constructs interaction terms
#' @param x matrix of predictors
#' @param top_ints k x 2 integer matrix of interaction terms
#' @seealso \code{\link[intscreen]{filter_top}}, \code{\link[intscreen]{compute_cors}}, and \code{\link[intscreen]{intscreen}}
#' @export
#' @examples
#'
#' x <- matrix(rnorm(1000 * 500), ncol = 500)
#' y <- rnorm(1000)
#'
#' cormat <- compute_cors(x, y)
#'
#' int_idx <- filter_top(cormat, k = 50)
#'
#' int_mat <- construct_ints(x, ints)
construct_ints <- function(x, top_ints)
{
    cnames <- colnames(x)

    if (is.null(cnames))
    {
        cnames <- paste0("X", 1:ncol(x))
    }

    intnames <- paste(cnames[top_ints[,1]], cnames[top_ints[,2]], sep = ":")

    int_mat <- construct_ints_cpp(X = x, whichints = top_ints)

    colnames(int_mat) <- intnames

    int_mat
}



#' all computations for interaction screening
#' @description This function does all computations for interaction screening
#' @param x matrix of predictors
#' @param y vector of observations
#' @param k integer number of top interactions to screen
#' @export
#' @examples
#'
#' x <- matrix(rnorm(1000 * 500), ncol = 500)
#' y <- rnorm(1000)
#'
#' ints <- intscreen(x, y, k = 12)
#'
#' str(ints)
intscreen <- function(x, y, k = 10)
{
    cormat <- compute_cors(x, y)

    int_idx <- filter_top(cormat, k = k)

    int_mat <- construct_ints(x, int_idx)

    list(int_mat = int_mat, int_idx = int_idx)
}


#' all computations for interaction screening
#' @description This function does all computations for interaction screening
#' @param x matrix of predictors
#' @param y vector of observations
#' @param k integer number of top interactions to screen
#' @export
#' @examples
#'
#' set.seed(1)
#' x <- matrix(rnorm(100 * 350), ncol = 350)
#' y <- rnorm(100) + x[,1] * x[,2] * 0.5 - x[,3] * x[,4] * 0.5
#'
#' ## require that each interaction be in the top 200 ints 95% of the 15 splits
#' ints <- cv_intscreen(x, y, k = 200, nsplits = 15, fraction.in.thresh = 0.95)
#'
#' ints$int_idx
#'
#' ## require that each interaction be in the top 50 ints 100% of the 15 splits
#' ints <- cv_intscreen(x, y, k = 50, nsplits = 15, fraction.in.thresh = 1)
#'
#' ints$int_idx
cv_intscreen <- function(x, y, k = 100, nsplits = 10, train.frac = 0.75, fraction.in.thresh = 1)
{
    intlist <- vector(mode = "list", length = nsplits)

    n <- nrow(x)

    for (s in 1:nsplits)
    {
        s_idx   <- sample.int(n, floor(n * train.frac))
        cormat  <- compute_cors(x[s_idx,,drop=FALSE], y[s_idx])
        intlist[[s]] <- filter_top(cormat, k = k)
    }

    ## the first index (ie first column of matrix returned by filter_top)
    ## is always less than the second index, so we don't need to worry
    ## about ordering of interactions
    int_tall <- do.call(rbind, intlist)

    all_ints <- paste(int_tall[,1], int_tall[,2], sep = ":")

    tab <- table(all_ints)

    ints_all_in <- tab[tab >= fraction.in.thresh * nsplits ]

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
        int_mat <- construct_ints(x, int_idx_final)
    }

    list(int_mat = int_mat, int_idx = int_idx_final)
}


