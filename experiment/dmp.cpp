#include "dmp.h"
#include <cmath>
#include <armadillo>
#include "scaling.cpp"


extern double dt;

/******** public ********/
dmp::dmp(uint nBF, double T, bool can_clock_index)
{
    this->clock = canclock(T,can_clock_index);
    this->w.zeros();
    this->regress = regressor(nBF);
}

arma::vec dmp::init(arma::mat y_data, arma::vec y0, arma::vec g)
{
    this->n_dof = y_data.n_rows;
    this->w = arma::zeros(n_dof,this->get_regressor().get_nBF());
    arma::vec error = this->trainDMP(y_data);
    
    this->y0d = y_data.col(0);
    this->gd = y_data.col(y_data.n_cols-1);
    this->y0 = y0;
    this->g = g;
    
    if ( this->n_dof == 3 )
        this->ksflag = true;
    else
        this->ksflag = false;

    this->setks();
    
    return error;
}

arma::vec dmp::getYx(arma::vec t)
{
    arma::vec x = this->clock.getPhase(t);
    arma::mat phi = this->regress.get_Fx(x);
    arma::vec yx;
    yx = (this->ks)*(this->w.t()*phi - this->y0d) + this->y0;
    return yx;
}

arma::vec dmp::getYx_dot(arma::vec t)
{
    arma::vec x = this->clock.getPhase(t);
    arma::vec x_dot = this->clock.getPhaseDot(t);
    arma::mat phi_dot = this->regress.get_dFx(x, x_dot);
    arma::vec yx_dot;
    yx_dot = (this->ks)*this->w.t()*phi_dot;
    return yx_dot;
}

arma::vec dmp::getYx_ddot(arma::vec t)
{
    arma::vec x = this->clock.getPhase(t);
    arma::vec x_dot = this->clock.getPhaseDot(t);
    arma::vec x_ddot = this->clock.getPhaseDDot(t);
    arma::mat phi_ddot = this->regress.get_ddFx(x, x_dot, x_ddot);
    arma::vec yx_ddot;
    yx_ddot = (this->ks)*this->w.t()*phi_ddot;
    return yx_ddot;
}

arma::mat dmp::simulation(double t, double dt, arma::vec y, arma::vec y_dot, arma::vec frep, double T)
{
    double K = 400;
    double D = 40;
    arma::vec yx = this->getYx(t*arma::ones(1));
    arma::vec yx_dot = this->getYx_dot(t*arma::ones(1));
    arma::vec yx_ddot = this->getYx_ddot(t*arma::ones(1));
    arma::vec ddy;
    ddy = yx_ddot - D*( y_dot - yx_dot ) - K*( y - yx ) + frep;
    arma::vec yr = y + y_dot*dt;
    arma::vec dyr = y_dot + ddy*dt;
    arma::mat Y = arma::zeros(3,3);
    Y.col(0) = yr;
    Y.col(1) = dyr;
    Y.col(2) = ddy;
    return Y;
}

void dmp::set_scaleMethod(bool flag)
{
    if (flag)
    {
        this->ksflag = 1;
        this->setks();
    }
    else
    {
        this->ksflag = 0;
        this->setks();
    }
}

// Setters.
void dmp::set_w(arma::mat w)
{
    this->w = w;
}

void dmp::set_g(arma::vec g)
{
    this->g = g;
}

void dmp::set_y0(arma::vec y0)
{
    this->y0 = y0;
}

void dmp::set_clock(canclock clock)
{
    this->clock = clock;
}

void dmp::setks()
{
    if ((!this->ksflag) || (this->n_dof != 3))
        this->ks = this->get1Dscaling();
    else
        this->ks = this->get3Dscaling();
}

//Getters.
arma::vec dmp::get_g()
{
    return this->g;
}

arma::vec dmp::get_gd()
{
    return this->gd;
}

arma::vec dmp::get_y0()
{
    return this->y0;
}

arma::vec dmp::get_y0d()
{
    return this->y0d;
}

regressor dmp::get_regressor()
{
    return this->regress;
}

canclock dmp::get_clock()
{
    return this->clock;
}

uint dmp::get_ndof()
{
    return this->n_dof;
}

arma::mat dmp::getW()
{
    return this->w;
}

arma::mat dmp::getks()
{
    return this->ks;
}

/******** private ********/

arma::vec dmp::trainDMP(arma::mat y_data)
{
    uint n = y_data.n_cols;
    arma::vec C = arma::linspace(0,1,n);
    arma::mat psi = this->regress.get_Fx(C);
    arma::vec error = this->set_weights(psi,y_data);
    return error;
}

arma::mat dmp::get1Dscaling()
{
    arma::mat ks = arma::zeros(3,3);
    ks.diag() = (this->g - this->y0)/(this->gd - this->y0d);
    return ks;
}

arma::mat dmp::get3Dscaling()
{
    return scaling(this->y0d, this->gd, this->y0, this->g);
}

arma::vec dmp::set_weights(arma::mat psi, arma::mat y_demo)
{
    this->w = arma::solve(psi.t(),y_demo.t());
    arma::mat error_mat = arma::pow(y_demo - arma::trans(this->w)*psi,2);
    arma::vec mse = arma::mean(error_mat,1);
    return mse;
}