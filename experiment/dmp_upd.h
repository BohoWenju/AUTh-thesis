#ifndef DMP_UPD_H
#define DMP_UPD_H

#include <armadillo>
#include "dmp.h"
#include "constraint.h"

class dmp_upd:public dmp
{
public:
    constraint con1; // Constraint to be removed in each update.
    constraint con2; // Constraint to parse in each update;
    arma::mat P; // Matrix for weight update.
    bool method; // Flag to showcase whether an optimization would happen on acceleration or position.
        // (Acelleration 1. Position 0.)
    bool adapt_online; // Flag. If it is 1 then the dmp will take the previous state
        // of the robot. For 0 it will take the refference signal as
        // previous state.
    double K; // Stiffness value.
    double D; // Damping value;
public:

    /* DMP_upd constructor 
    ~  @param[in] nBF: The number of basis functions(kernels).
    ~  @param[in] T : Scale for canonical clock.
    ~  @param[in] can_clock_index : Index to choose which canonical clock to use[0:exp,1:linear].   
    ~  @param[in] method : Flag to choose optimization. (Acelleration 0. Position 1.).
    ~ For the rest of dmp_upd parameters simply the null constructor of the constraint
    ~ class should be called. The __init handles the initialization.
    */
    dmp_upd(uint nBF, double T, bool can_clock_index, bool method, constraint con2);
    dmp_upd(uint nBF, double T);


    /* Initialization function. Handles the training. And the first update.
        ~  @param[in] y_data: Demonstrated position trajectory.
        ~  @param[in] y0 : Desired initial position.
        ~  @param[in] g : Desired final position(goal). 
        ~  @param[out] error : Returns training error.
    */
    arma::vec init_upd(arma::mat y_data, double dt, arma::vec y0, arma::vec g);
    arma::vec init_upd(arma::mat y_data, double dt);

    /*  Function that simulates the system.
        ~  @param[in] t : Current time moment.
        ~  @param[in] y : Previous position.
        ~  @param[in] y_dot : Previous velocity.
        ~  @param[in] online : Online flag. Set to 1 when there are changes.
        ~  @param[in] con : New constraint.
        ~  @param[in] frep : External force implemented as coupling term.
        ~  @param[in] T : In case the simulation time differs from the one
        ~  in training data.
        ~  @param[in] a_f : In case the scaling of the yx_ddot term changes
        ~  online.a E [0,1] : For a = 0 the system behaves normally.
        ~  For a E (0,1] the system adapts better to previous states.
        ~  @param[in] camera_flag : When an updated position from camera
        ~  comes.
        ~  @param[out] Y : Matrix in the form of : Y = [y dy ddy].
        */
    arma::mat simulation(double t, double dt, arma::vec y, arma::vec y_dot, bool online,
    constraint con, arma::vec frep, double T, bool camera_flag);

    arma::mat simulation(double t, double dt, arma::vec y, arma::vec y_dot, bool online,
    constraint con, arma::vec frep, double T);

    arma::mat simulation(double t,double dt, arma::vec y, arma::vec y_dot, bool online,
    constraint con, arma::vec frep);

    constraint get_con2();

public:
    /* Function that calculates the scaling of acceleration.
        ~ @param[in] frep : coupling term.
    */
    double calc_af(arma::vec frep);
    /* Initialization for constraints.
       The necessary values are in the dmp super class.
    */
    void initConstraints();

    // Initialize P matrix.
    void initP();

    /* Function that transforms position based on on h_transform.
        ~ @param[in] y : Position to apply h_transform.
        ~ @param[in] flag : Flag to choose between scaling method.
    */
    arma::vec h_transform(arma::vec y, bool flag);

    /* Function to change the scaling term.
        ~ @param[in] con : Constraint from which the new goal and initial position will be extracted.
    */
    void update_ks(constraint con);

    /* Function to update the constraints.
        ~ @param[in] con : New Constraint.
    */
    void updConstraints(constraint con);

    /* Update. This function updates the weights of the dmp by removing the previous
     constraints and altering the weights so as to fulfil the current constraints.
     This is done by solving the KKT problem. The code below is derived from the eq. (HERE)
        ~ @param[in] t : Current moment. This is needed as one necessary constraint is the
      remainder of the previous state. If t == -1 then the initial update occurs.
        ~ @param[in] con : Desired new constraint.
        ~ @param[in] camera_flag : Flag to augment the importance of
      updating goal constraint.
    */
    void update(double t, double dt, constraint con, bool camera_flag);

    /*Weight Update. This function updates the weights of the dmp and the P matrix.
        ~  @param[in] e : Sensitivity matrix. Holds the sensitivity for each constraint.
        ~  @param[in] flag : Flag that showcases whether to remove the constraint or add it.
    */
    void update_weights(arma::mat Z, arma::mat H, arma::vec e, bool flag);
       
    // The functions below are used for the weight upd. No explanation is needed.
    arma::mat get_Phat(arma::mat K, arma::mat H, arma::mat R);

    arma::mat get_What(arma::mat K, arma::mat Z, arma::mat H);

    arma::mat get_Khat(arma::mat H, arma::mat R);

    /* These functions returns matrices Z,H and also sensitivity matrix e.
      As input the current constraint is parsed in order to extract from it the 
      necessary constraints. Also, time moment t is needed to form the necessary constraints.
        ~  @param[in] t : Time moment.
        ~  @param[in] dt : Needed to extract previous acceleration.
        ~  @param[in] con : Constraint in order to extract necessary values.
        ~  @param[in] flag : Flag to know whether an update or a downdate occurs.
                        1 : update, 0 : downdate.
        ~  @param[in] camera_flag : Explained above.
    */
    arma::vec get_e(double t, constraint con, bool flag, bool camera_flag);

    arma::mat get_H(double t, double dt, constraint con);

    arma::mat get_Z(double t, double dt, constraint con);

    /* Get A matrix.
       This function returns A matrix which is a matrix that contains basis function values
       for a given moment t. Also it containts the time derivatives of the basis functions.
        ~  @param[in] t : Time moment.
    */
   arma::mat get_A(double t);

    /* Function to transform current constraints based on h_transform.
     Also inverse transform is aplied.
        ~  @param[in] flag : Regular or inverse h_transform
    */
   void con_transform(bool flag);

   /* Function to get reference position for dmp. When optimization in acceleration
     occurs in order to get the position one only needs: yr = W'*phi. No scaling is needed.
        ~  @param[in] t : Time moment.
    */
   arma::vec getYr(double t);

    /* Velocity of reference signal.
        ~ @param[in] t : Time moment.
    */
   arma::vec getYr_dot(double t);

    /* Acceleration of reference signal.
        ~ @param[in] t : Time moment.
    */
   arma::vec getYr_ddot(double t);
};
#endif