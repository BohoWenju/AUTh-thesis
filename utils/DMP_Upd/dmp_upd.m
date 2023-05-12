classdef dmp_upd < dmp
    methods (Access = public)
        % DMP_upd constructor 
        % For the dmp constructor:
        %  @param[in] nBF: The number of basis functions(kernels).
        %  @param[in] T : Scale for canonical clock.
        %  @param[in] can_clock_index : Index to choose which canonical clock to use[0:exp,1:linear].   
        %  @param[in] method : Flag to choose optimization. (Acelleration 0. Position 1.).
        % For the rest of dmp_upd parameters simply the null constructor of the constraint
        % class should be called. The __init handles the initialization.
        function this = dmp_upd(nBF,T,can_clock_index,method)
            this@dmp(nBF,T,can_clock_index);
            if (nargin < 4)
                this.method = 0;
            else
                this.method = method;
            end
            this.con1 = constraint();
            this.con2 = constraint();
            this.set_adapt(0);
        end

        % Initialization function
        %  @param[in] T : Timesteps that correspond to y_data.
        %  @param[in] y_data : Demonstrated position trajectory.
        %  @param[in] y0 : Desired initial position
        %  @param[in] g : Desired final position(goal).
        %  The above parameters also hold the initial desired constraints.
        function init_upd(this,T,y_data,y0,g)
            % Calling the init of the superclass
            this.init(T,y_data,y0,g);
            % Initializing constraints.
            this.initConstraints();
            this.initP();
            % Initial update on weights.
            if this.method
                % Apply h_transform to constraints.
                this.con_transform(1);
                % Weights need to be transformed as well. 
                this.set_w(this.getW()*(this.getks()'));
            end
            
            [Z, H, e] = this.getZHe(-1,this.con2);

            % Update based on new constraints.
            this.update_weights(Z,H,e,1);

            if this.method
                % Apply inverse h_transform to constraints.
                this.con_transform(0);
                % Inverse transform for weights.
                this.set_w(this.getW()/(this.getks()'));
            end
        end

        % Function that simulates the system.
        %  @param[in] t : Current time moment.
        %  @param[in] y : Previous position.
        %  @param[in] y_dot : Previous velocity.
        %  @param[in] online : Online flag. Set to 1 when there are changes.
        %  @param[in] con : New constraint.
        %  @param[in] frep : External force implemented as coupling term.
        %  @param[in] T : In case the simulation time differs from the one
        %  in training data.
        %  @param[in] a_f : In case the scaling of the yx_ddot term changes
        %  online.a E [0,1] : For a = 0 the system behaves normally.
        %  For a E (0,1] the system adapts better to previous states.
        %  @param[in] camera_flag : When an updated position from camera
        %  comes.
        function [y,dy,ddy] = simulation(this,t,y,y_dot,online,con,frep,T,a_f,camera_flag)
            
            % Global step of integration
            global dt
            c_flag = 0;
            a = 0.1;
            if ( nargin < 7 )
                frep = 0;
                a = 0.1;
            elseif ( nargin < 9 )
                this.set_clock(canClock(T));
                a = 0.1;
            elseif ( nargin < 10 )
                a = a_f;
                this.set_clock(canClock(T));              
            else
                a = a_f;
                this.set_clock(canClock(T));
                c_flag = camera_flag;
            end
            
            if online
                if (~this.adapt_online)
                    yx = this.getYx(t);
                    yx_dot = this.getYx_dot(t);
                    if t > 0
                        yx_ddot = this.getYx_ddot(t-dt);
                    else
                        yx_ddot = this.getYx_ddot(0);
                    end
                    Y = [yx yx_dot yx_ddot];
                    con.set_ypr(Y);
                end
                % Update Weights
                this.update(t,con,c_flag);
                this.update_ks(con);
            end
            
            % Desired stiffness parameter
            K = 300;
            % Desired damping parameter
            D = 2*(sqrt(K+10));

            if this.method
                % Get reference values
                yx = this.getYx(t);
                yx_dot = this.getYx_dot(t);
                yx_ddot = this.getYx_ddot(t);
            else
                % Get reference values
                yx = this.getYr(t);
                yx_dot = this.getYr_dot(t);
                yx_ddot = this.getYr_ddot(t);
            end

            % Integration
            ddy = (1-a)*yx_ddot - D*(y_dot-yx_dot) - K*(y-yx) + frep;
            y = y + y_dot*dt;
            dy = y_dot+ddy*dt;
        end
        
        % Function to set adapt online flag.
        function set_adapt(this,flag)
            this.adapt_online = flag;
        end
    end
    

    
    methods (Access = private)
        % Initialization for constraints.
        % The necessary values are in the dmp super class.
        function initConstraints(this)
%             this.con1.con_upd(this.get_gd);
            this.con2.con_upd(this.get_g);
            this.con2.set_y0(this.get_y0);
        end

        % Initialize P matrix. Initialization is derived from eq. 
        % P matrix initialization does not need any particular input.
        % (HERE)
        function initP(this)
            x_data = 0:0.01:1;
            x_dot = 1;
            x_ddot = 0;
            n_data = length(x_data);
            H = zeros(this.get_regressor.get_nBF, n_data);
            if this.method
                for j=1:n_data,H(:,j) = this.get_regressor.get_Fx(x_data(j));end
            else
                for j=1:n_data,H(:,j) = this.get_regressor.get_ddFx(x_data(j),x_dot,x_ddot);end
            end
            tol = .1;
            this.P = inv(H*H' + tol*eye(this.get_regressor.get_nBF));
        end

        % Function that transforms position based on on h_transform.
        % The code below is based on eq. (HERE)
        % Also, the inverse h_tranform is available based on a flag parsed as input.
        %  @param[in] y : Position to apply h_transform. If flag = 0 then y is the output
        %  of the h_transform.
        %  @param[in] flag : Flag to choose between h_transform or inverse h_transform.
        function h = h_transform(this,y,flag)
            if flag
                h = (y-this.get_y0())+this.getks()*this.get_y0d();
            else
                h = y + this.get_y0() -this.getks()*this.get_y0d();
            end
        end

        % Function to update scaling term due to target or initial position
        % change.
        %  @param[in] con : New constraint from which we want to extract
        %  the new position and goal.
        function update_ks(this,con)
            if con.isGoal
                this.set_g(con.get_goal());
                if con.isY0
                    this.set_y0(con.get_y0());
                end
                this.setks();
            else
                if con.isY0
                    this.set_y0(con.get_y0());
                    this.setks();
                end
            end
        end
                
        % Update constraints. This function has different structure from the init_constraints
        % because when a new constraint is added the previous constraint that was supposed
        % to be fulfilled must be the one to be removed.
        %  @param[in] con : Desired constraint.
        function updConstraints(this,con)
            this.con1.copy(this.con2);
            this.con1.isYpr = 0;
            this.con1.isY0 = 0;
            this.con2.copy(con);
        end

        
        % Update. This function updates the weights of the dmp by removing the previous
        % constraints and altering the weights so as to fulfil the current constraints.
        % This is done by solving the KKT problem. The code below is derived from the eq. (HERE)
        %  @param[in] t : Current moment. This is needed as one necessary constraint is the
        %  remainder of the previous state. If t == -1 then the initial update occurs.
        %  @param[in] con : Desired new constraint.
        %  @param[in] camera_flag : Flag to augment the importance of
        %  updating goal constraint.
        function update(this,t,con,camera_flag)
            
            this.updConstraints(con);
            
            if this.method
                % Apply h_transform to constraints.
                this.con_transform(1);
                % Weights need to be transformed as well. 
                this.set_w(this.getW()*(this.getks()'));
            end

            % Downdate occurs only at runtime
            % Remove previous constraints
            [Z , H, e] = this.getZHe(t,this.con1,0,camera_flag);
            this.update_weights(Z,H,e,0);
            
            % Get new constraints.
            [Z, H, e] = this.getZHe(t,this.con2,1,camera_flag);

            % Update based on new constraints.
            this.update_weights(Z,H,e,1);

            if this.method
                % Apply inverse h_transform to constraints.
                this.con_transform(0);
                % Inverse transform for weights.
                this.set_w(this.getW()/(this.getks()'));
            end
        end


        % Weight Update. This function updates the weights of the dmp and the P matrix.
        % Most of the input parameters defined here are explained in eq. (HERE)
        %  @param[in] e : Sensitivity matrix. Holds the sensitivity for each constraint.
        %  @param[in] flag : Flag that showcases whether to remove the constraint or add it.
        function update_weights(this,Z,H,e,flag)
            if (e < 0)
                return;
            end
            if ~flag
                e = -e;
            end
            R = diag(e);
            K = get_Khat(this,H,R);
            this.set_w(this.get_What(K,Z,H));
            this.P = this.get_Phat(K,H,R);
        end

        % The functions below are used for the weight upd. No explanation is needed.
        % eq. (HERE)
        function P = get_Phat(this,K,H,R)
            I = eye(size(this.P))-H*K;
            P = I'*this.P*I+K'*R*K;
            P = (P+P')/2;
        end
        
        % eq. (HERE)
        function W = get_What(this,K,Z,H)
            W = this.getW()+((Z-(this.getW()')*H)*K)';
        end
        
        % eq. (HERE)
        function K = get_Khat(this,H,R)
%             K = (this.P*H)/(R+(H')*this.P*H);
            C = this.P*H;
            dS = decomposition(H'*this.P*H + R, 'ldl', 'upper');
            K = (dS \ C');
        end


        % This function returns matrices Z,H and also sensitivity matrix e.
        % As input the current constraint is parsed in order to extract from it the 
        % necessary constraints. Also, time moment t is needed to form the necessary constraints.
        %  @param[in] t : Time moment.
        %  @param[in] con : Constraint in order to extract necessary values.
        %  @param[in] flag : Flag to know whether an update or a downdate occurs.
        %                    1 : update, 0 : downdate.
        %  @param[in] camera_flag : Explained above.
        function [Z,H,e] = getZHe(this,t,con,flag,camera_flag)
    
            persistent A0 A1;
            n = this.get_ndof;
            Z = zeros(n,1);
            H = zeros(this.get_regressor.get_nBF,1);
            change_flag = 0;
            n_constraints = 1;
            e(:,1) = zeros(3,1);
            
            if t < 0   
                A0 = this.get_A(0);
                A1 = this.get_A(this.get_clock.get_tau);
                if con.isYpr 
                    con.set_y0(con.get_ypr);
                    con.isYpr = 0;
                end
            end
            
            % Check if previous states are concerned.
            if con.isYpr
                [Z,H] = this.setPr(t,con,Z,H);
                e_t = 1e3*[1e-6; 1e-4; 1e-3];
                if (length(e_t))~=1
                    for i=1:length(e_t(1,:))
                        e(:,i) = e_t(:,i);
                        n_constraints = n_constraints+1;
                    end
                end
                change_flag=1;
            end
            
            
            % Check if change in initial Position is added.
            if con.isY0
                if t < 0 
                    Y0 = [con.get_y0 zeros(n,1) zeros(n,1)];
                    Z = [Z Y0];
                    H = [H A0];
                    e(:,n_constraints) = [1e-9; 1e-7; 1e-7];
                else
                    Y0 = con.get_y0;
                    Z = [Z Y0];
                    H = [H A0(:,1)];
                    if flag
                        e(1,n_constraints) = 1e-9;
                    else
                        e(1,n_constraints) = 1e-9;
                    end
                end
                n_constraints = n_constraints +1;
                change_flag = 1;
            end
            
            
            % Check if change in final Position is added.
            if con.isGoal
                if t < 0 
                    Y0 = [con.get_goal zeros(n,1) zeros(n,1)];
                    Z = [Z Y0];
                    H = [H A1];
                    e(:,n_constraints) =  1e6*[1e-1; 1e-0; 1e-0];
                else
                    Y0 = con.get_goal;
                    Z = [Z Y0 zeros(3,1) zeros(3,1)];
                    H = [H A1];
                    if flag
                        if camera_flag
                            e(:,n_constraints) = 1e0*[1e-9; 1e-7; 1e-7];
                        else
                            e(:,n_constraints) = 1e6*[1e-1; 1e-0; 1e-0];
                        end
                    else
                        if camera_flag
                            e(:,n_constraints) = 1e6*[1e-1; 1e-0; 1e-0];
                        else
                            e(:,n_constraints) = 1e6*[1e-1; 1e-0; 1e-0];
                        end
                    end
                end   
                n_constraints = n_constraints +1;
                change_flag =1;
            end
         
            % Check if via point is added.
            if con.isVia
                via = con.get_via;
                t_via = cons.t_via;
                x_via = this.get_clock.getPhase(t_via);
                psi = this.get_regressor.get_Fx(x_via);
                Z = [Z via(:,i)];
                H = [H psi];
                if flag
                    e(:,n_constraints) = 1e-7;
                else
                    e(:,n_constraints) = 1e-7;
                end
                n_constraints = n_constraints +1;
                change_flag = 1;
            end
            
            
            if n_constraints ~= 1
                e = matrix2vec(e);
            else
                e = -1;
            end
                
            if change_flag
                Z = Z(:,2:end);
                H = H(:,2:end);
            end
        end

        % This function sets Z and H matrix according to [psi(t) dpsi(t) ddpsi(t-dt)].
        % The parameter con should contain the corresponding values.
        %  @param[in] t : Time moment.
        %  @param[in] con : Constraint in order to extract previous values.
        %  @param[in] Z,H : Matrices that will be augmented.
        function [Z,H] = setPr(this,t,con,Z,H)
            persistent A0;
            if t > 0
                A = this.get_A(t);
                Ci = [A(:,1:2) A0(:,3)];
                A0 = Ci;
            else
                Ci = this.get_A(0);
                A0 = Ci;
            end 
            Z = [Z con.get_ypr];
            H = [H Ci];
        end

        % Get A matrix.
        % This function returns A matrix which is a matrix that contains basis function values
        % for a given moment t. Also it containts the time derivatives of the basis functions.
        %  @param[in] t : Time moment.
        function A = get_A(this,t)
            x = this.get_clock.getPhase(t);
            x_dot = this.get_clock.getPhaseDot(t);
            x_ddot = this.get_clock.getPhaseDDot(t);
            psi = this.get_regressor.get_Fx(x);
            dpsi = this.get_regressor.get_dFx(x,x_dot);
            ddpsi = this.get_regressor.get_ddFx(x,x_dot,x_ddot);
            A = [psi dpsi ddpsi];
        end

        % Function to transform current constraints based on h_transform.
        % Also inverse transform is aplied.
        %  @param[in] flag : Regular or inverse h_transform.
        function con_transform(this,flag)
            % For first constrain.
            if this.con1.isGoal
                this.con1.set_goal(this.h_transform(this.con1.get_goal(),flag));
            end
            if this.con1.isY0
                this.con1.set_y0(this.h_transform(this.con1.get_y0(),flag));
            end
            if this.con1.isYpr
                this.con1.set_ypr(this.h_transform(this.con1.get_ypr(),flag));
            end
            if this.con1.isVia
                this.con1.set_ypr(this.h_transform(this.con1.get_via(),flag));
            end

            % For second constraint.
            if this.con2.isGoal
                this.con2.set_goal(this.h_transform(this.con2.get_goal(),flag));
            end
            if this.con2.isY0
                this.con2.set_y0(this.h_transform(this.con2.get_y0(),flag));
            end
            if this.con2.isYpr
                this.con2.set_ypr(this.h_transform(this.con2.get_ypr(),flag));
            end
            if this.con2.isVia
                this.con2.set_ypr(this.h_transform(this.con2.get_via(),flag));
            end


        end

        % Function to get reference position for dmp. When optimization in acceleration
        % occurs in order to get the position one only needs: yr = W'*phi. No scaling is needed.
        %  @param[in] t : Time moment.
        function yr = getYr(this,t)
            x = this.get_clock.getPhase(t);
            yr = (this.getW)'*this.get_regressor.get_Fx(x);
        end

        % Velocity of reference signal.
        %  @param[in] t : Time moment.
        function yr = getYr_dot(this,t)
            x = this.get_clock.getPhase(t);
            x_dot = this.get_clock.getPhaseDot(t);
            yr = (this.getW)'*this.get_regressor.get_dFx(x,x_dot);
        end

        % Acceleration of reference signal.
        %  @param[in] t : Time moment.
        function yr = getYr_ddot(this,t)
            x = this.get_clock.getPhase(t);
            x_dot = this.get_clock.getPhaseDot(t);
            x_ddot = this.get_clock.getPhaseDDot(t);
            yr = (this.getW)'*this.get_regressor.get_ddFx(x,x_dot,x_ddot);
        end

    end

    properties (Access = private)
        % This object inherits the properties of the super class dmp.
        con1 % Constraints to be removed in each update.
        con2 % Constraints to parse in each update.
        P % P matrix for weight update. 
        method % Whether optimization on acceleration lever or position will occur.
        % (Acelleration 1. Position 0.).
        adapt_online % Flag. If it is 1 then the dmp will take the previous state
        % of the robot. For 0 it will take the refference signal as
        % previous state.
    end
end

