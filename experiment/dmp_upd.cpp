#include "dmp_upd.h"

#include <cmath>
#include <armadillo>


/******** public ********/
dmp_upd::dmp_upd(uint nBF, double T, bool can_clock_index, bool method, constraint con2):dmp(nBF,T,can_clock_index)
{
    this->method = method;
    this->con1 = constraint();
    this->con2 = con2;
    this->adapt_online = true;
    this->K = 300;
    this->D = 2*std::sqrt(K);
}

dmp_upd::dmp_upd(uint nBF, double T):dmp(nBF,T,true)
{
    this->method = 0;
    this->con1 = constraint();
    this->con2 = constraint();
    this->adapt_online = true;
    this->K = 300;
    this->D = 2*std::sqrt(K);
}

arma::vec dmp_upd::init_upd(arma::mat y_data, double dt, arma::vec y0, arma::vec g)
{
    dmp::set_g(g);
    dmp::set_y0(y0);
    arma::vec error = dmp::init(y_data,y0,g);
    this->initConstraints();
    this->initP();
    arma::mat Z = this->get_Z(-1,dt,this->con2);
    arma::mat H = this->get_H(-1,dt,this->con2);
    arma::vec e = this->get_e(-1,this->con2,true,false);
    
    this->update_weights(Z,H,e,1);
    return error;
}

arma::vec dmp_upd::init_upd(arma::mat y_data, double dt)
{
    dmp::set_g(y_data.col(y_data.n_cols-1));
    dmp::set_y0(y_data.col(0));
    arma::vec error = dmp::init(y_data,dmp::get_y0(),dmp::get_g());
    
    this->initConstraints();
    
    this->initP();

    arma::mat Z = this->get_Z(-1,dt,this->con2);
    arma::mat H = this->get_H(-1,dt,this->con2);
    arma::vec e = this->get_e(-1,this->con2,true,false);
    
    this->update_weights(Z,H,e,1);
    return error;
}

double dmp_upd::calc_af(arma::vec frep)
{
    double fmax = 1e-3;
    double temp = arma::norm(frep)/fmax;
    arma::vec temp1 = {temp,1};
    return  temp1.min();
}

arma::mat dmp_upd::simulation(double t, double dt, arma::vec y, arma::vec y_dot, bool online,
constraint con, arma::vec frep, double T, bool camera_flag)
{
    double a_f = this->calc_af(frep);
    arma::vec yx = this->getYr(t);
    arma::vec yx_dot = this->getYr_dot(t);
    arma::vec yx_ddot = this->getYr_ddot(t);
    arma::mat Y;
    dmp::set_clock(canclock(T,dmp::get_clock().get_index()));
    if (online) 
    {
        if (this->adapt_online)
        {
            arma::vec yx_ddotCon;
            if (t > 0)
                yx_ddotCon = this->getYr_ddot(t-dt);
            else
                yx_ddotCon = this->getYr_ddot(0);
            Y = arma::join_horiz(yx,arma::join_horiz(yx_dot,yx_ddotCon));
            con.set_ypr(Y);
        }
        this->update(t,dt,con,camera_flag);
    }

    
    arma::vec ddy = (1-a_f)*yx_ddot - this->D*(y_dot-yx_dot) - this->K*(y-yx) + frep;
    y = y + y_dot*dt;
    arma::vec dy = y_dot+ddy*dt;
    arma::mat Y1 = arma::join_horiz(y,join_horiz(dy,ddy));
    return Y1;
}

arma::mat dmp_upd::simulation(double t, double dt, arma::vec y, arma::vec y_dot, bool online,
constraint con, arma::vec frep, double T)
{
    double a_f = this->calc_af(frep);
    arma::vec yx = this->getYr(t);
    arma::vec yx_dot = this->getYr_dot(t);
    arma::vec yx_ddot = this->getYr_ddot(t);
    arma::mat Y;
    dmp::set_clock(canclock(T,dmp::get_clock().get_index()));
    if (online) 
    {
        if (this->adapt_online)
        {
            arma::vec yx_ddotCon;
            if (t > 0)
                yx_ddotCon = this->getYr_ddot(t-dt);
            else
                yx_ddotCon = this->getYr_ddot(0);
            Y = arma::join_horiz(yx,arma::join_horiz(yx_dot,yx_ddotCon));
            con.set_ypr(Y);
        }
        this->update(t,dt,con,false);
    }

    arma::vec ddy = (1-a_f)*yx_ddot - this->D*(y_dot-yx_dot) - this->K*(y-yx) + frep;
    y = y + y_dot*dt;
    arma::vec dy = y_dot+ddy*dt;
    arma::mat Y1 = arma::join_horiz(y,join_horiz(dy,ddy));
    return Y1;
}

arma::mat dmp_upd::simulation(double t, double dt, arma::vec y, arma::vec y_dot, bool online,
constraint con, arma::vec frep)
{
    if (online) 
    {
        arma::mat Y;
        if (this->adapt_online)
        {
            arma::vec yx = this->getYr(t);
            arma::vec yx_dot = this->getYr_dot(t);
            arma::vec yx_ddot = this->getYr_ddot(t);
            arma::vec yx_ddotCon;
            if (t > 0)
                yx_ddotCon = this->getYr_ddot(t-dt);
            else
                yx_ddotCon = this->getYr_ddot(0);
            Y = arma::join_horiz(yx,arma::join_horiz(yx_dot,yx_ddotCon));
            con.set_ypr(Y);
        }
        this->update(t,dt,con,false);
        return Y;
    }

    double a_f = this->calc_af(frep);
    arma::vec yx = this->getYr(t);
    arma::vec yx_dot = this->getYr_dot(t);
    arma::vec yx_ddot = this->getYr_ddot(t);

    arma::vec ddy = (1-a_f)*yx_ddot - this->D*(y_dot-yx_dot) - this->K*(y-yx) + frep;
    y = y + y_dot*dt;
    arma::vec dy = y_dot+ddy*dt;
    arma::mat Y1 = arma::join_horiz(y,join_horiz(dy,ddy));
    return Y1;
}


constraint dmp_upd::get_con2()
{
    return this->con2;
}

/******** private ********/
void dmp_upd::initConstraints()
{
    this->con2.con_upd(dmp::get_g());
    this->con2.set_y0(dmp::get_y0());
}

void dmp_upd::initP()
{
    uint n = 100;
    arma::vec x_data = arma::linspace(0,1,n);
    arma::mat H = arma::zeros(dmp::get_regressor().get_nBF(),n);
    for (int j=0; j<n; j++)
        H.col(j) = dmp::get_regressor().get_ddFx(x_data(j)*arma::ones(1,1),arma::ones(1,1),arma::zeros(1,1));
    double tol = 0.1;
    arma::mat temp;
    temp.eye(dmp::get_regressor().get_nBF(),dmp::get_regressor().get_nBF());
    this->P = H*H.t() + tol*temp;
    this->P = solve(this->P, temp);
}

arma::vec dmp_upd::h_transform(arma::vec y, bool flag)
{
    arma::vec h;
    if (flag)
        h = (y-dmp::get_y0()) + dmp::getks()*dmp::get_y0d();
    else
        h = y + dmp::get_y0() - dmp::getks()*dmp::get_y0d();
    return h;
}

void dmp_upd::update_ks(constraint con)
{
    if (con.isGoal)
        dmp::set_g(con.get_goal());

    if (con.isY0)
        dmp::set_y0(con.get_y0());
    dmp::setks();
}

void dmp_upd::updConstraints(constraint con)
{
    this->con1.copy(this->con2);
    this->con1.isYpr = false;
    this->con1.isY0 = false;
    this->con2.copy(con);
}

void dmp_upd::update(double t, double dt, constraint con, bool camera_flag)
{
    this->updConstraints(con);
    arma::mat Z = this->get_Z(t,dt,this->con1);
    arma::mat H = this->get_H(t,dt,this->con1);
    arma::vec e = this->get_e(t,this->con1,true,false);

    this->update_weights(Z,H,e,0);

    arma::mat Z1 = this->get_Z(t,dt,this->con2);
    arma::mat H1 = this->get_H(t,dt,this->con2);
    arma::vec e1 = this->get_e(t,this->con2,true,false);
    this->update_weights(Z1,H1,e1,1);
}

void dmp_upd::update_weights(arma::mat Z, arma::mat H, arma::vec e, bool flag)
{
    arma::mat R = arma::diagmat(e);
    arma::mat K = this->get_Khat(H,R);
    dmp::set_w(this->get_What(K,Z,H));
    this->P = this->get_Phat(K,H,R);
}

arma::mat dmp_upd::get_Phat(arma::mat K, arma::mat H, arma::mat R)
{
    arma::mat I = arma::eye(this->P.n_rows, this->P.n_cols);
    arma::mat P = arma::zeros(this->P.n_rows, this->P.n_cols);
    I = I - H*K;
    P = I.t()*this->P*I+K.t()*R*K;
    P = (P+P.t())/2;
    return P;
}

arma::mat dmp_upd::get_What(arma::mat K, arma::mat Z, arma::mat H)
{
    arma::mat W;
    W = dmp::getW() + arma::trans((Z-(arma::trans(dmp::getW()))*H)*K);
    return W;
}

arma::mat dmp_upd::get_Khat(arma::mat H, arma::mat R)
{
    arma::mat C = this->P*H;
    arma::mat K = arma::inv(R+H.t()*this->P*H);
    return arma::solve(K,C.t());
}

arma::vec dmp_upd::get_e(double t, constraint con, bool flag, bool camera_flag)
{
    arma::mat e = arma::zeros(1,1);
    bool change_flag = false;
    if (t < 0)
        con.isYpr = false;

    if (con.isYpr)
    {
        arma::vec e_t = {1e-4,1e-1,1e-1};
        e  = arma::join_vert(e,e_t);
        change_flag = true;
    }

    if (con.isY0)
    {
        arma::vec e_t;
        if (t < 0)
            e_t = {1e-9,1e-7,1e-7};
        else
            if (flag) // update
            {
                if (camera_flag)
                    e_t = {1e-4,1e-1,1e-1};
                else
                    e_t = {1e+3, 1e+9, 1e+9};
            }
            else // downdate
                e_t = {1e+3, 1e+9, 1e+9};
        e  = arma::join_vert(e,e_t);
        change_flag = true;
    }

    if (con.isGoal)
    {
        arma::vec e_t;
        if (t < 0)
            e_t = {1e-9,1e-7,1e-7};
        else
        if (flag) // update
        {
            if (camera_flag)
                e_t = {1e-4,1e-1,1e-1};
            else
                e_t = {1e+3, 1e+9, 1e+9};
        }
        else // downdate
            e_t = {1e+3, 1e+9, 1e+9};
        e  = arma::join_vert(e,e_t);
        change_flag = true;
    }
    if (change_flag)
        e = e.rows(1,e.n_rows-1);
    return e;
}

arma::mat dmp_upd::get_H(double t, double dt, constraint con)
{
    uint n = dmp::get_ndof();
    arma::mat H = arma::zeros(dmp::get_regressor().get_nBF(),1); 
    bool change_flag = false;
    if (t < 0)
        con.isYpr = false;

    if (con.isYpr)
    {
        arma::mat A, A1;
        if (t > 0)
        {
            A = this->get_A(t);
            A1 = this->get_A(t-dt);
        }
        else
        {
            A = this->get_A(0);
            A1 = this->get_A(0);
        }
        H = arma::join_horiz(H,A.cols(0,1),A1.col(2));
        change_flag = true;
    }
    if (con.isY0)
    {
        arma::mat A0 = this->get_A(0);
        H = arma::join_horiz(H, A0);
        change_flag = true;
    }

    if (con.isGoal)
    {
        arma::mat A0 = this->get_A(dmp::get_clock().get_tau());
        H = arma::join_horiz(H, A0);
        change_flag = true;
    }

    // no via point implementation since it is not needed.
    if (change_flag) 
        H = H.cols(1,H.n_cols-1);
    return H;
}

arma::mat dmp_upd::get_Z(double t, double dt, constraint con)
{
    uint n = dmp::get_ndof();
    arma::mat Z = arma::zeros(n,1); 
    uint n_cosntraints = 0;
    bool change_flag = false;
    if (t < 0)
        con.isYpr = false;

    if (con.isYpr)
    {
        Z = arma::join_horiz(Z,con.get_ypr());
        change_flag = true;
    }
    if (con.isY0)
    {
        Z = arma::join_horiz(Z,con.get_y0());
        Z = arma::join_horiz(Z,arma::zeros(Z.n_rows,2));
        change_flag = true;
    }

    if (con.isGoal)
    {
        Z = arma::join_horiz(Z,con.get_goal());
        Z = arma::join_horiz(Z,arma::zeros(Z.n_rows,2));
        change_flag = true;
    }

    // no via point implementation since it is not needed.
    if (change_flag){
        uint n_cols = Z.n_cols;
        Z = Z.cols(1,n_cols-1);
    }
    return Z;
}

arma::mat dmp_upd::get_A(double t)
{
    arma::vec x = (dmp::get_clock()).getPhase(t*arma::ones(1,1));
    arma::vec x_dot = (dmp::get_clock()).getPhaseDot(t*arma::ones(1,1));
    arma::vec x_ddot = (dmp::get_clock()).getPhaseDDot(t*arma::ones(1,1));
    arma::vec psi = (dmp::get_regressor()).get_Fx(x);
    arma::vec dpsi = (dmp::get_regressor()).get_dFx(x,x_dot);
    arma::vec ddpsi = (dmp::get_regressor()).get_ddFx(x,x_dot,x_ddot);
    arma::mat A = arma::ones(dmp::get_regressor().get_nBF(),3);
    A.col(0) = psi;
    A.col(1) = dpsi;
    A.col(2) = ddpsi;
    return A;
}

arma::vec dmp_upd::getYr(double t)
{
    arma::vec x = (dmp::get_clock()).getPhase(t*arma::ones(1,1));
    arma::vec yr = arma::trans((dmp::getW()))*(dmp::get_regressor()).get_Fx(x);
    return yr;
}

arma::vec dmp_upd::getYr_dot(double t)
{
    arma::vec x = (dmp::get_clock()).getPhase(t*arma::ones(1,1));
    arma::vec x_dot = (dmp::get_clock()).getPhaseDot(t*arma::ones(1,1));
    arma::vec dyr = arma::trans((dmp::getW()))*(dmp::get_regressor()).get_dFx(x,x_dot);
    return dyr;
}

arma::vec dmp_upd::getYr_ddot(double t)
{
    arma::vec x = (dmp::get_clock()).getPhase(t*arma::ones(1,1));
    arma::vec x_dot = (dmp::get_clock()).getPhaseDot(t*arma::ones(1,1));
    arma::vec x_ddot = (dmp::get_clock()).getPhaseDDot(t*arma::ones(1,1));
    arma::vec ddyr = arma::trans((dmp::getW()))*(dmp::get_regressor()).get_ddFx(x,x_dot,x_ddot);
    return ddyr;
}

