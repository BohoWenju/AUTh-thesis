#include <armadillo>
arma::mat get_skew(arma::vec a)
{
    arma::mat temp = arma::zeros(3,3);
    temp(0,1) = -a(2);
    temp(0,2) = a(1);
    temp(1,0) = a(2);
    temp(1,2) = -a(0);
    temp(2,0) = -a(1);
    temp(2,1) = a(0);
    return temp;
}