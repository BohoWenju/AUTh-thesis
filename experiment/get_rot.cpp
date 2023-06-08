#include <armadillo>
#include <cmath>
#include "get_skew.cpp"
arma::mat get_rot(arma::vec n, arma::vec nd)
{
    double cost = arma::sum(n%nd); 
    arma::vec k = get_skew(nd)*n;
    double sint = arma::norm(k);
    arma::mat Rg;
    arma::mat Sg = get_skew(k);
    Rg = arma::eye(3,3) + Sg*sint + Sg*Sg*(1-cost);
    return Rg;
}