#include "controller.h"
#include <armadillo>
#include <cmath>

controller::controller(arma::mat Kpos, arma::mat Krot)
{
    this->Kpos = Kpos;
    this->Krot = Krot;
}

arma::mat controller::simulation(double t, double dt, arma::mat xe, arma::vec yd,
     arma::vec yd_dot, arma::vec q, arma::vec dq, arma::vec fext)
{

arma::mat Je = jacobian(q);
arma::mat Re = this.robot.fkine(q);
arma::mat Re_temp = Re(1:3,1:3);
Re_temp = Re_temp.t();
arma::vec eo = this.pdRot_control(q,yd_dot);

% Closed loop inverse kinematics.
arma::vec up = yd_dot - this.Kpos*(xe - yd);
arma::vec ur = +this.Krot*eo;
arma::vec u = arma::joint_horiz(up,ur);
arma::vec dq1 = arma::solve(Je,u);
% Integration.
arma::vec q = q + dq1*dt;
% Derivation.
arma::vec d2q = (dq1-dq)/dt;
arma::mat return_matrix = arma::joint_horiz(q,arma::joint_horiz(dq,d2q));
}

arma::vec controller::pdRot_control(arma::vec q, arma::vec dy)
{
    arma::mat temp = this.robot.fkine(q);
    arma::vec x_axis = make_unitary(temp(1:3,1));
    arma::vec ndy = make_unitary(dy);
    if (ismequal(ndy,zeros(3,1)))
        ndy = n;
    end
    arma::vec ne = get_skew(x_axis)*ndy;
    double theta = atan2(norm(ne),x_axis.t()*ndy);
    return theta*ne;
}