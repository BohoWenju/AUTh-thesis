#ifndef CANCLOCK_H
#define CANCLOCK_H

#include <armadillo>

class canclock
{
private:
    bool index; // Flag to choose between desired canonical clock.
    double tau; // Scaling parameter.
    double ac; // Parameter for exponential canonical clock.
public:
    /* Constructor for class canclock.
        ~ @param[in] scale : Parameter to define the scaling of the canonical system.
        ~ @param[in] can_clock_index : Flag to define whether the exponential canonical clock(0) should be used or the linear (~0)
    */
    canclock(double scale, bool can_clock_index = 1);

    // Null contructor.
    canclock();
    
    /* Returns the phase based on a current moment t.
        ~ @param[in] t : Moment(s) on which the desired phase value should be calculated.
    */
    arma::vec getPhase(arma::vec t);

    /* Returns the time derivative of current phase t.
        ~ @param[in] t : Moment(s) on which the desired phase_dot value should be calculated.
    */
    arma::vec getPhaseDot(arma::vec t);

    /* Returns the second time derivative of current phase t.
        ~ @param[in] t : Moment(s) on which the desired phase_ddot value should be calculated.
    */
    arma::vec getPhaseDDot(arma::vec t);

    // Getter for scale parameter.
    double get_tau();

    bool get_index();
};

#endif
