#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <armadillo.h>
class controller
{
private:
    arma::mat Kpos;
    arma::mat Krot:
public:
    controller(arma::mat Kpos, arma::mat Krot);
    arma::mat simulation(double t, double dt, arma::mat xe, arma::vec yd,
     arma::vec yd_dot, arma::vec q, arma::vec qd, arma::vec fext);
    arma::vec pd_rotControl(arma::vec q, arma::vec dy);
};



#endif 