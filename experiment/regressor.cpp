#include "regressor.h"

#include <cmath>
#include <armadillo>

/******** public ********/
regressor::regressor(double nBF)
{
    this->nBF = nBF;
    this->set_centers();
    this->set_variance(1.5);
}

regressor::regressor()
{
    this->nBF = 40;
    this->set_centers();
    this->set_variance(1.5);
}

arma::mat regressor::get_Fx(arma::vec x)
{
    uint16_t n = x.n_rows;
    arma::mat psi = this->get_psi(x);
    arma::mat Fx = arma::mat(this->nBF,n);
    double spsi;
    for (int i = 0; i < n; i++)
    {
        if ( x(i) < 0 ) {psi.col(i) = arma::ones(this->nBF,i);}
        else if ( x(i) > 1 ) {psi.col(i) = arma::zeros(this->nBF,i);}
        else 
        {
            spsi = arma::sum(psi.col(i));
            Fx.col(i) = psi.col(i)/spsi;
        }
    }
    return Fx;
}

arma::mat regressor::get_dFx(arma::vec x, arma::vec x_dot)
{
    uint16_t n = x.n_rows;
    arma::mat psi = this->get_psi(x);
    arma::mat dpsi = this->get_psiDot(x,x_dot);
    arma::mat Fx = this->get_Fx(x);
    arma::mat dFx = arma::mat(this->nBF,n);
    double spsi;
    double dspsi;
    for (int i = 0; i < n; i++)
    {
        if ( x(i) < 0 ) {psi.col(i) = arma::zeros(this->nBF,i);}
        else if ( x(i) > 1 ) {psi.col(i) = arma::zeros(this->nBF,i);}
        else 
        {
            spsi = arma::sum(psi.col(i));
            dspsi = arma::sum(dpsi.col(i));
            dFx.col(i) = (dpsi.col(i)-Fx.col(i)*dspsi)/spsi;
        }
    }    
    return dFx;
}

arma::mat regressor::get_ddFx(arma::vec x, arma::vec x_dot, arma::vec x_ddot)
{
    uint16_t n = x.n_rows;
    arma::mat psi = this->get_psi(x);
    arma::mat dpsi = this->get_psiDot(x,x_dot);
    arma::mat Fx = this->get_Fx(x);
    arma::mat dFx = this->get_dFx(x,x_dot);
    arma::mat ddpsi = this->get_psiDDot(x,x_dot,x_ddot);
    arma::mat ddFx = arma::mat(this->nBF,n);
    double spsi;
    double dspsi;
    double ddspsi;
    for (int i = 0; i < n; i++)
    {
        if ( x(i) < 0 ) {psi.col(i) = arma::zeros(this->nBF,i);}
        else if ( x(i) > 1 ) {psi.col(i) = arma::zeros(this->nBF,i);}
        else 
        {
            spsi = arma::sum(psi.col(i));
            dspsi = arma::sum(dpsi.col(i));
            ddFx.col(i) = (ddpsi.col(i) - 2*dFx.col(i)*dspsi - Fx.col(i)*ddspsi)/spsi;
        }
    }
    return ddFx;
}

arma::vec regressor::get_c()
{
    return this->c;
}

double regressor::get_h()
{
    return this->h;
}

uint16_t regressor::get_nBF()
{
    return this->nBF;
}

/******** protected ********/
arma::mat regressor::get_psi(arma::vec x)
{
    uint16_t n = x.n_rows;
    arma::mat psi = arma::mat(this->nBF,n);
    for (int i = 0; i < n; i++)
        psi.col(i) = arma::exp(-this->h*(arma::pow(x(i)-this->c,2)));

    return psi;
}

arma::mat regressor::get_psiDot(arma::vec x, arma::vec x_dot)
{
    uint16_t n = x.n_rows;
    arma::mat psi = arma::mat(this->nBF,n);
    psi = this->get_psi(x);
    arma::mat dpsi = arma::mat(this->nBF,n);
    arma::mat a = arma::mat(this->nBF,n);
    for (int i = 0; i < n; i++)
        a.col(i) = ( x(i) - this->c)*x_dot(i);

    dpsi = -2*this->h*(psi % a);
    return dpsi;
}

arma::mat regressor::get_psiDDot(arma::vec x, arma::vec x_dot, arma::vec x_ddot)
{
    uint16_t n = x.n_rows;
    arma::mat psi = this->get_psi(x);
    arma::mat dpsi = this->get_psiDot(x,x_dot);
    arma::mat ddpsi = arma::mat(this->nBF,n);
    arma::mat a = arma::mat(this->nBF,n);
    arma::mat da = arma::mat(this->nBF,n);
    for (int i =0; i < n; i++)
    {
        a.col(i) = ( x(i) - this->c)*x_dot(i);
        da.col(i) = (x(i) - this->c)*x_ddot(i) + std::pow(x_dot(i),2);
    }

    ddpsi = -2*this->h*((dpsi % a) + (psi % da));
    return ddpsi;
}

/******** private ********/
void regressor::set_centers()
{
    this->c = arma::linspace(0,1,this->nBF);
}

void regressor::set_variance(double ah)
{
    double temp = this->c(1) - this->c(0);
    temp = std::pow((ah*temp),2);
    this->h = 1/temp;
}

