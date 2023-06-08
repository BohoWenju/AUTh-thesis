#include "canclock.h"

#include <cmath>
#include <armadillo>

/******** public ********/
canclock::canclock(double scale, bool can_clock_index)
{
    this->index = can_clock_index;
    this->tau = scale;
    this->ac = 1;
}

canclock::canclock()
{
    this->index = 1;
    this->tau = 10;
    this->ac = 1;
}

arma::vec canclock::getPhase(arma::vec t)
{
    arma::vec x = arma::vec(t.n_rows);
    if (this->index){x = (1/this->tau)*t;}
    else {x = exp(-(t*this->ac)/this->tau);}
    return x;
}

arma::vec canclock::getPhaseDot(arma::vec t)
{
    arma::vec dx = arma::vec(t.n_rows);
    if (this->index){dx = (1/this->tau)*arma::ones(t.n_rows);}
    else {dx = -(this->ac/this->tau)*exp(-(t*this->ac)/this->tau);}
    return dx;
}

arma::vec canclock::getPhaseDDot(arma::vec t)
{
    arma::vec ddx = arma::vec(t.n_rows);
    if (this->index){ddx.zeros();}
    else {ddx = +(pow(this->ac,2)/pow(this->tau,2))*exp(-(t*this->ac)/this->tau);}
    return ddx;
}

double canclock::get_tau()
{
    return this->tau;
}

bool canclock::get_index()
{
    return this->index;
}


