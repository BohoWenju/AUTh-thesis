#include <iostream>
#include <armadillo>
#include "dmp_upd.h"
using namespace std;
using namespace arma;
double dt = 0.002;

cube get5thTrajectory(vec y0, vec yf, vec Time)
{
    double T = Time(Time.n_rows-1);
    vec t = Time/T;
    mat y;
    y = (yf - y0)*(trans(10*arma::pow(t,3) - 15*arma::pow(t,4) + 6*arma::pow(t,5)));
    for(int i = 0; i<Time.n_rows; i++)
      y.col(i) += y0;
    mat y_dot;
    y_dot = (yf - y0)*(trans(30*arma::pow(t,2) - 60*arma::pow(t,3) + 30*arma::pow(t,4)))/T;
    mat y_ddot;
    y_ddot = (yf - y0)*(trans(60*t - 180*arma::pow(t,2) + 120*arma::pow(t,3)))/(T*T);
    cube C(y.n_rows,y.n_cols,3);
    C.slice(0) = y;
    C.slice(1) = y_dot;
    C.slice(2) = y_ddot;
    return C;
}

int main()
{

  double dt = 0.002;
  // Initialize training data.
  vec y0 = zeros(3,1);
  vec yf = ones(3,1);
  cube C = get5thTrajectory(y0,yf,linspace(0,10,10));
  mat y_data = C.slice(0);

  // Initialize constraint.
  constraint con = constraint(yf,zeros(3,3),y0, zeros(3,1));
  
  dmp_upd temp = dmp_upd(10,10);
  arma::mat A = temp.init_upd(y_data,dt);
  
  mat Y = ones(3,3);
  con.set_ypr(Y);
  con.set_goal(2*yf);
  con.isVia = 0;
  con.isY0 = 0;
  arma::mat temp_sim = temp.simulation(0,dt,y0,zeros(3,1),true,con,zeros(3,1));
  temp.updConstraints(con);
  cout << temp_sim <<endl;
  return 0;
}


