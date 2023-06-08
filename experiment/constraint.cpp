#include "constraint.h"

#include <armadillo>

/******** public ********/
constraint::constraint(arma::vec g, arma::mat ypr, arma::vec y0,arma::vec via)
{
    this->set_goal(g);
    this->set_ypr(ypr);
    this->set_y0(y0);
    this->set_via(via);
}

constraint::constraint(){
    
}

void constraint::con_upd(arma::vec g, arma::mat ypr, arma::vec y0,arma::vec via)
{
    this->set_goal(g);
    this->set_ypr(ypr);
    this->set_y0(y0);
    this->set_via(via);
}

void constraint::con_upd(arma::vec g, arma::mat ypr, arma::vec y0)
{
    this->set_goal(g);
    this->set_ypr(ypr);
    this->set_y0(y0);
    this->isVia = false;
}

void constraint::con_upd(arma::vec g, arma::mat ypr)
{
    this->set_goal(g);
    this->set_ypr(ypr);
    this->isY0 = false;
    this->isVia = false;
}

void constraint::con_upd(arma::vec g)
{
    this->set_goal(g);
    this->isY0 = false;
    this->isYpr = false;
    this->isVia = false;
}
void constraint::deactivate(bool goalflag, bool y0flag, bool yprflag, bool viaflag)
{
    this->isGoal = goalflag;
    this->isY0 = y0flag;
    this->isYpr = yprflag;
    this->isVia = viaflag;
}

void constraint::copy(constraint con)
{
    if (con.isGoal)
        this->set_goal(con.get_goal());
    else
        this->isGoal = 0;
    if (con.isY0)
        this->set_y0(con.get_y0());
    else
        this->isY0 = 0;
    if (con.isYpr)
        this->set_ypr(con.get_ypr());
    else
        this->isYpr = 0;
    if (con.isVia)
        this->set_via(con.get_via());
    else
        this->isVia = 0;
}

void constraint::remove()
{
    this->isYpr = 0;
}

arma::vec constraint::get_goal()
{
    arma::vec g;
    if (this->isGoal)
        g = this->g;
    else
        throw std::runtime_error("Goal is not parsed");
    return g;
}

arma::vec constraint::get_y0()
{
    arma::vec y0;
    if (this->isY0)
        y0 = this->y0;
    else
        throw std::runtime_error("Y0 is not parsed");
    return y0;
}

arma::mat constraint::get_ypr()
{
    arma::mat ypr;
    if (this->isYpr)
        ypr = this->ypr;
    else
        throw std::runtime_error("Ypr is not parsed");
    return ypr;
}

arma::mat constraint::get_via()
{
    arma::mat via;
    if (this->isVia)
        via = this->via;
    else
        throw std::runtime_error("Via is not parsed");
    return via;
}

void constraint::set_goal(arma::vec g)
{
    this->g = g;
    this->isGoal = 1;
}

void constraint::set_y0(arma::vec y0)
{
    this->y0 = y0;
    this->isY0 = 1;
}

void constraint::set_ypr(arma::mat ypr)
{
    this->ypr = ypr;
    this->isYpr = 1;
}

void constraint::set_via(arma::mat via)
{
    this->via = via;
    this->isVia = 1;
}