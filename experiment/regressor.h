#ifndef REGRESSOR_H
#define REGRESSOR_H

#include <armadillo>

class regressor
{
private:
    double nBF; // Number of kernels.
    arma::vec c; // Centers for kernels.
    double h; // Variance of kernel functions.

    /* Function that sets the centers of the kernels. Centers should be equally spaced in time between 0 ~ 1.
    No parameter is parsed in this function due to the nBF parameter already in class.
    */
    void set_centers();

    /* Function that sets the variance of the kernels. Must be called after set_centers.
    No parameter is parsed in this function due to the c parameter already in class.
    */
    void set_variance(double ah = 1.5);

protected:
    /* Get Psi values. Numerator of Fx components.
        ~ @param[in] x : Canonical clock output for a moment t. It is a vector due to the possibility
        of multiple time moments.
    */
   arma::mat get_psi(arma::vec x);

   /* Get dPsi/dt values. d(Numerator of Fx components)/dt.
        ~ @param[in] x : Canonical clock output for a moment t. It is a vector due to the possibility
        of multiple time moments.
        ~ @param[in] x_dot : Canonical clock time derivative output.
    */
   arma::mat get_psiDot(arma::vec x, arma::vec x_dot);

   /* Get ddPsi/ddt values. dd(Numerator of Fx components)/ddt.
        ~ @param[in] x : Canonical clock output for a moment t. It is a vector due to the possibility
        of multiple time moments.
        ~ @param[in] x_dot : Canonical clock time derivative output.
         ~ @param[in] x_ddot : Canonical clock second time derivative output.
    */
   arma::mat get_psiDDot(arma::vec x, arma::vec x_dot, arma::vec x_ddot);
     
public:
    /* Constructor for class regressor.
        ~ @param[in] nBF : Number of kernels.
    */
    regressor(double nBF);

    // Null constructor.
    regressor();
    
    /* Get Fx values.
        ~ @param[in] x : Canonical clock output for a moment t. It is a vector due to the possibility
        of multiple time moments.
    */
    arma::mat get_Fx(arma::vec x);

    /* Get dFx/dt values.
        ~ @param[in] x : Canonical clock output for a moment t. It is a vector due to the possibility
        of multiple time moments.
        ~ @param[in] x_dot : Canonical clock time derivative output.
    */
    arma::mat get_dFx(arma::vec x, arma::vec x_dot);

    /* Get ddFx/dt values.
        ~ @param[in] x : Canonical clock output for a moment t. It is a vector due to the possibility
        of multiple time moments.
        ~ @param[in] x_dot : Canonical clock time derivative output.
        ~ @param[in] x_ddot : Canonical clock second time derivative output.
    */
    arma::mat get_ddFx(arma::vec x, arma::vec x_dot, arma::vec x_ddot);

    // Getter for centers.
    arma::vec get_c();

    // Getter for variance.
    double get_h();

    // Getter for nBF.
    uint16_t get_nBF();
};

#endif