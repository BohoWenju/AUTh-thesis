#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <armadillo>

class constraint
{
private:
    arma::vec g; // Desired goal.
    arma::vec y0; // Desired initial position.
    arma::mat ypr; // Previous state of the DMP. Necessary for online updates.
    arma::mat via; // Desired via points.
public:
    bool isGoal; // Flag that states if the goal is an active constraint.
    bool isY0; // Flag that states if the inital position is an active constraint.
    bool isYpr; // Flag that states if the previous state is an active constraint.
    bool isVia; // Flag that states if the via point is an active constraint.
    
    /* Constructor for class Constraint.
        ~  @param[in] g : Desired target.
        ~  @param[in] y0 : Desired initial position.
        ~  @param[in] ypr : Previous state of the DMP.
        ~  @param[in] via : Desired via point(s).
    */
    constraint(arma::vec g, arma::mat ypr, arma::vec y0, arma::vec via);

    //Null constructor 
    constraint();
    /* Function that updates the constraints based on new desired target,initial position
        etc.
       This function has the same input and structure as the constructor 
        for obvious reasons.
    */
    void con_upd(arma::vec g, arma::mat ypr, arma::vec y0, arma::vec via);
    void con_upd(arma::vec g, arma::mat ypr, arma::vec y0);
    void con_upd(arma::vec g, arma::mat ypr);
    void con_upd(arma::vec g);


    /* Function that sets the desired flags to desired values.
        ~  @param[in] goalflag : Goal flag.
        ~  @param[in] y0flag : Initial position flag.
        ~  @param[in] yprflag : Previous state flag.
        ~  @param[in] viaflag : Via point(s) flag.
    */
    void deactivate(bool goalflag, bool y0flag, bool yprflag, bool viaflag);

    /* Function that copies another constraint.
        ~  @param[in] con : constraint to copy.
    */
   void copy(constraint con);

   /* Function that removes Ypr.*/
   void remove();

   // Getters.
   arma::vec get_goal();

   arma::mat get_ypr();

   arma::mat get_via();

   arma::vec get_y0();

   // Setters. Obvious parameters no need for explanation.
   void set_goal(arma::vec goal);

   void set_ypr(arma::mat ypr);

   void set_y0(arma::vec y0);

   void set_via(arma::mat via);

};

#endif
