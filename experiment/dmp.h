#ifndef DMP_H
#define DMP_H

#include <armadillo>
#include "canclock.h"
#include "regressor.h"
class dmp
{
private:
    arma::mat w; // Weights of the DMP.
    arma::mat ks; // Scaling term.
    uint n_dof; // Degree of freedom.
    canclock clock; // Canonical clock.
    regressor regress; // Regressor object to obtain regreessor values.
    arma::vec g; // Desired goal.
    arma::vec y0; // Desired initial position.
    arma::vec gd; // Goal extracted from demo.
    arma::vec y0d; // Initial position extracted from demo.
    bool ksflag; // Flag to determine scaling parameter.
public:

    /* DMP Constructor.
        ~  @param[in] nBF: The number of basis functions(kernels).
        ~  @param[in] T : Scale for canonical clock.
        ~  @param[in] can_clock_index : Index to choose which canonical clock to use[0:exp,1:linear]. 
    */
    dmp(uint nBF, double T, bool can_clock_index);

    /* DMP Initialization function.
        ~  @param[in] y_data: Demonstrated position trajectory.
        ~  @param[in] y0 : Desired initial position.
        ~  @param[in] g : Desired final position(goal). 
        ~  @param[out] error : Returns training error.
    */
    arma::vec init(arma::mat y_data, arma::vec y0, arma::vec g);

    /* Function that returns the reference position trajectory.
        ~  @param[in] t : Current time moment. Can be a number of moments in the case of a vector.
    */ 
    arma::vec getYx(arma::vec t);

    /* Function that returns the reference position trajectory.
        ~  @param[in] t : Current time moment. Can be a number of moments in the case of a vector.
    */ 
    arma::vec getYx_dot(arma::vec t);

    /* Function that returns the reference position trajectory.
        ~  @param[in] t : Current time moment. Can be a number of moments in the case of a vector.
    */ 
    arma::vec getYx_ddot(arma::vec t);

    /* Function that simulates the system.
        ~  @param[in] t : Current time moment. 
        ~  @param[in] dt : Timestep of integration. 
        ~  @param[in] y : Current position.
        ~  @param[in] y_dot : Current velocity.
        ~  @param[in] frep : External force applied as a coupling term.
        ~  @param[in] T : In case simulation time changes then clock scale should also change.
        ~  @param[out] Y : Matrix in the form of : Y = [y dy ddy].
    */ 
    arma::mat simulation(double t, double dt, arma::vec y, arma::vec y_dot, arma::vec frep, double T);

    /* Function that sets the desired scaling method.
        ~  @param[in] flag : Desired scale method.
    */
    void set_scaleMethod(bool flag);

    // Getters.
    arma::vec get_g();

    arma::vec get_gd();

    arma::vec get_y0();

    arma::vec get_y0d();

    regressor get_regressor();

    canclock get_clock();

    uint get_ndof();

    arma::mat getW();

    arma::mat getks();

    // Setters.
    void set_w(arma::mat w);

    void set_g(arma::vec g);

    void set_y0(arma::vec y0);

    void set_clock(canclock clock);

    void setks();

private:
    /* Function that trains the DMP weights on a trajectory.
        ~  @param[in] y_data : Demo position trajectory.
        ~  @param[out] error : Training error.
    */
    arma::vec trainDMP(arma::mat y_data);

    // Function that sets the scaling term based on eq11 from this paper : A Reversible Dynamic Movement Primitive formulation.
    arma::mat get1Dscaling();

    // Function that sets the scaling term based on eq11 from this paper : A novel DMP formulation for global and frame independent spatial scaling in the task space.
    arma::mat get3Dscaling();

    /* Function that trains the DMP weights on a trajectory.
        ~  @param[in] psi : Regressor values.
        ~  @param[in] y_demo : Demonstrated trajectory.
        ~  @param[out] mse : Mean squared error.
    */
    arma::vec set_weights(arma::mat psi, arma::mat y_demo);
};





#endif