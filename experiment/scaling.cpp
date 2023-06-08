#include <cmath>
#include <armadillo>
#include "get_rot.cpp"
arma::mat scaling(arma::vec y0d, arma::vec gd, arma::vec y0, arma::vec g)
{
    double temp = norm(g-y0);
    double temp1 = norm(gd-y0d);
    arma::vec nd = (gd-y0d)/temp1;
    arma::vec n = (g-y0)/temp;

    
    double sg = temp/temp1;
    arma::mat Rg = get_rot(n, nd);
    arma::mat ks = sg*Rg;
    return ks;
}
