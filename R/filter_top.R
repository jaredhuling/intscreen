
#' filters top interactions using correlation matrix
#' @description This function filters top interactions using correlation matrix
#' @param cormat correlation matrix returned by \code{compute_cors}
#' @param k integer number of top interactions to screen
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
#' @export
compute_cors <- function(x, y)
{
    compute_cors_cpp(X = x, Y = y)
}


#' constructs selected interactions
#' @description This function constructs interaction terms
#' @param x matrix of predictors
#' @param top_ints k x 2 integer matrix of interaction terms
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
#' cormat <- compute_cors(x, y)
#'
#' int_idx <- filter_top(cormat, k = 50)
#'
#' int_mat <- construct_ints(x, y, k = 12)
#'
#' head(int_mat)
screen_ints <- function(x, y, k = 10)
{
    cormat <- compute_cors(x, y) / NROW(x)

    int_idx <- filter_top(cormat, k = k)

    int_mat <- construct_ints(x, int_idx)

    list(int_mat = int_mat, int_idx = int_idx)
}

