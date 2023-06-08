#include <iostream>
#include <armadillo>
#include "dmp.h"

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

  // Initialize training data.
  vec y0 = zeros(3,1);
  vec yf = ones(3,1);
  cube C = get5thTrajectory(y0,yf,linspace(0,10,10));
  mat y_data = C.slice(0);
  
  // Initialize dmp object. Train object.
  dmp temp = dmp(40,10,true);
  vec error = temp.init(y_data,y0,yf);

  //Simulation parameters.
  double T = 10;
  double dt = 0.002;
  cube Y = zeros(3,3,((int)(T/dt))+1);
  int i = 0;
  vec y = y0;
  vec y_dot = zeros(3,1);
  mat exit_y = zeros(3,Y.n_slices);
  mat temp_mat = zeros(3,3);
  // Simulation.
  
  for (double t = 0; t < T; t+=dt)
  {
    temp_mat = temp.simulation(t,dt,y,y_dot,zeros(3,1),T);
    cout<<temp_mat.col(0)<<endl;
    y = temp_mat.col(0);
    y_dot = temp_mat.col(1);
    i++;
  }
  
  return 0;
}


