#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

//[[Rcpp::export]]
Eigen::MatrixXd compute_cors_cpp(const Eigen::Map<Eigen::MatrixXd> & X, const Eigen::Map<Eigen::VectorXd> & Y)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    int p = X.cols();

    MatrixXd cormat(MatrixXd::Zero(p, p));

    for (int i = 0; i < p - 1; ++i)
    {
        for (int j = i + 1; j < p; ++j)
        {
            cormat(i,j) = Y.dot((X.col(i).array() * X.col(j).array()).matrix() );
        }
    }

    return cormat;
}

//[[Rcpp::export]]
Eigen::MatrixXd construct_ints_cpp(const Eigen::Map<Eigen::MatrixXd> & X, const Eigen::Map<Eigen::MatrixXi> & whichints)
{
    using Eigen::MatrixXd;

    int p = X.cols();
    int n = X.rows();

    int num_ints = whichints.rows();

    MatrixXd intmat(n, num_ints);

    for (int i = 0; i < num_ints; ++i)
    {
        int var1 = whichints(i,0);
        int var2 = whichints(i,1);

        if (var1 > p || var2 > p)
        {
            stop("invalid variable index");
        }

        intmat.col(i) = X.col(var1).array() * X.col(var2).array();
    }

    return intmat;
}
