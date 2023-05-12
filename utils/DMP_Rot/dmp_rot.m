classdef dmp_rot < handle
    methods (Access = public)
        % Note : This class does not inherit dmp properties due to obvious reasons.
        % This class is based on : A correct formulation for the Orientation Dynamic Movement Primitives for robot control in the Cartesian space.
        % DMP rotation constructor.
        %  @param[in] nBF: The number of basis functions(kernels).
        %  @param[in] T : Scale for canonical clock.
        %  @param[in] az : With an appropriate choice of az, bz the system of the dmp becomes critically damped.
        %  @param[in] bz : The particular choice is : az = 4*bz
        %  @param[in] can_clock_index : Index to choose which canonical clock to use[0:exp,1:linear].
        function this = dmp_rot(nBF,T,can_clock_index,az,bz)
            if (nargin < 3)
                this.clock = canClock(T,1);
                this.az = 15;
                this.bz = this.az/4;
            elseif (nargin < 4)
                this.clock = canClock(T,can_clock_index);
                this.az = 15;
                this.bz = this.az/4;
            elseif (nargin < 5)
                this.clock = canClock(T,can_clock_index);
                this.az = az;
                this.bz = az/4;
            else
                this.clock = canClock(T,can_clock_index);
                this.az = az;
                this.bz = bz;
            end
            this.w = zeros(nBF,3);
            this.regressor = regressor(nBF);
            this.tau = T;
        end

        % Initialization function
        %  @param[in] Qd : Demonstrated quaternion trajectory.
        %  @param[in] T : Timesteps that correspond to orientation data.
        %  @param[in] w : Demonstrated angular velocity trajectory.
        function eqd = init(this,Qd,T,w)
            global dt
            this.Q0 = Qd(:,1);
            this.Qg = Qd(:,end);
            this.Q0d = Qd(:,1);
            this.Qgd = Qd(:,end);
            if (nargin < 3)
                T = this.clock.get_tau();
                T = 0:dt:T;
                [eqd, eqd_dot, eqd_ddot] = this.extract_error(Qd);
            elseif (nargin < 4)
                [eqd, eqd_dot, eqd_ddot] = this.extract_error(Qd);
            else
                [eqd, eqd_dot, eqd_ddot] = this.extract_error(Qd,w);
            end
            this.trainDMP(T,eqd, eqd_dot, eqd_ddot);
        end 

        % Simulation function.
        % @param[in] t : Current moment.
        % @param[in] eq : Orientation error.
        % @param[in] eq_dot : Time derivative of orientation error.
        % @param[in] w_dot : Angular acceleration.
        % @param[in] coupling_term : Coupling term in the form of a torque to be added in the eq. (14)
        % The final 2 parameters are for the extreme case of the goal or initial orientation changing 
        % during the runtime.
        %  @param[in] Qd : Demonstrated quaternion trajectory.
        %  @param[in] Q0 : Initial orientation.
        function [Q, w, w_dot, eq, eq_dot, eq_ddot] = simulation(this,t,eq,eq_dot,w,coupling_term,Qg,Q0)
            if (nargin < 6) %% PROBLEM HERE. WRITE IT BETTER
                coupling_term = zeros(3,1);
            elseif ( nargin < 7)  
            elseif ( nargin < 8)
                this.Qg = Qg;
            elseif (nargin < 9)
                this.Q0 = Q0;
                this.Qg = Qg;
            end
            persistent i
            if t == 0
                i = 1;
            end
            global dt
            % Formulating eq. (14)
            % needed for the forcing term.
            x = this.clock.getPhase(t);
            psi = this.regressor.get_Fx(x); 
%             this.set_scaling();
            eq_ddot = -this.az*(this.bz*eq+this.tau*eq_dot) + this.Fx(:,i) + coupling_term;
            eq_ddot = eq_ddot/(this.tau^2);
            eq = eq + eq_dot*dt;
            eq_dot = eq_dot + eq_ddot*dt;

            Q = prod_quat(bar_quat(exp_quat((1/2)*eq)),this.Qg);
            Q_dot = this.getQ_dot(Q,eq_dot);
            % Getting w from eq. (16)
            temp = prod_quat(Q_dot,bar_quat(Q));
            wn = 2*temp(2:4); % This is w(t+dt)
            w_dot = (wn - w)/dt;
            w = wn;
            i = i + 1;
        end


        % Function that sets the scaling term.
        function set_scaling(this)
            temp = 2*log_quat(prod_quat(this.Qg,bar_quat(this.Q0)));
            temp1 = 2*log_quat(prod_quat(this.Qgd,bar_quat(this.Q0d)));
            temp = diag(temp);
            temp1 = diag(temp1);
            this.S = temp1\temp;
        end
        
    end

    methods (Access = private)
        % Train the DMP
        %  @param[in] T : Timesteps that correspond to y_data.
        %  @param[in] eqd : Quaternion error between current orientation and goal orientation values.
        %  @param[in] eqd_dot : Time derivative of eqd.
        %  @param[in] eqd_ddot : Time derivative of eqd_dot.
        function trainDMP(this,T,eqd, eqd_dot, eqd_ddot)
            C = this.clock.getPhase(T);
            psi = this.regressor.get_Fx(C);

            % Calculate right side of eq. (23)
            this.set_Fx(eqd,eqd_dot,eqd_ddot);
            this.set_weights(psi.*C);
        end
        
        % Function that calculates the non-linear forcing term used for training the weights.
        %  @param[in] eqd : Quaternion error between current orientation and goal orientation values.
        %  @param[in] eqd_dot : Time derivative of eqd.
        %  @param[in] eqd_ddot : Time derivative of eqd_dot.
        function Fx = set_Fx(this,eqd,eqd_dot,eqd_ddot)
%             tempd = diag(2*log_quat(prod_quat(this.Qgd,bar_quat(this.Q0d))));
            tempd = eye(3);
            temp = (this.tau^2)*eqd_ddot + this.az*(this.bz*eqd + this.tau*eqd_dot);
            Fx = tempd\temp;
            this.Fx = Fx;
        end

        % Function that adapts the weights of the dmp kernels based on least squares objective function.
        % @param[in] psi : Regressor values. This vector should be nBF x m where m is the total number of timestamps.
        % @param[in] yd : Demonstrated trajectory.
        % @param[out] w : Weights for the trajectory derived from Least Squares.
        % The weights should be : nBF x n_dof (for formula correspondance)
        function set_weights(this,psi)
              yd = this.Fx;
            this.w = ((yd*(psi'))/(psi*(psi')))';
%               this.w = linsolve(psi',yd');
        end

        % Function that returns the quaternion error between the demonstrated trajectory and the goal quaternion.
        % Arithmetic derivation is used. This is purely symbolic due to the fact that the robot does not have velocity or acceleration sensors.
        % @param[in] Qd : Demonstrated orientation trajectory in the form of quaternions.
        % @param[in] w : Demonstrated angular velocity.
        function [eq, eq_dot, eq_ddot] = extract_error(this,Qd,w)
            n = length(Qd(1,:));
            for i = 1:n
                eq(:,i) = 2*log_quat(prod_quat(this.Qg,bar_quat(Qd(:,i))));
            end

            % Arithmetical derivative.
            global dt
            n = length(Qd(1,:));
            if (nargin < 3)
                for i = 1:n
                    if i + 1 > n
                        eq_dot(:,i) = eq_dot(:,i-1);
                        eq_ddot(:,i) = eq_ddot(:,i-1);
                        continue; % Exiting the loop.
                    elseif i + 2 > n
                        eq_dot(:,i) = (eq(:,i+1) - eq(:,i)) / dt;
                        eq_ddot(:,i) = eq_ddot(:,i-1);
                    else
                        eq_dot(:,i) = (eq(:,i+1) - eq(:,i)) / dt; 
                        eq_ddot(:,i) = (eq(:,i+2)-2*eq(:,i+1) + eq(:,i)) / (dt^2); 
                    end
                end
            else
                % Get Q_dot out of angular velocity.
                for i = 1:n
                    Q_dot = (1/2)*get_JQ(Qd(:,i))*w(:,i);
                    % Get eq_dot.
                    eq_dot(:,i) = this.get_deq(Qd(:,i),Q_dot);
                end
                
                % Arithmetical derivative for eq_ddot.
                for i = 1:n
                    if i + 1 > n
                        eq_ddot(:,i) = eq_ddot(:,i-1);
                    else
                        eq_ddot(:,i) = (eq_dot(:,i+1)-eq_dot(:,i))/dt;
                    end
                end
            end
        end

        % Function that returns the time derivative of a quaternion in order o calculate angular velocity.
        % Combination of eq. (21) and eq. (16)
        % @param[in] Q : Quaternion at current moment.
        % @param[in] eq_dot : Time derivative of orientation error.
        function Q_dot = getQ_dot(this,Q,eq_dot)
            temp = prod_quat(this.Qg,bar_quat(Q)); % Qg x Q_bar
            J = getJlog(temp);
            temp1 = J*eq_dot; % JlogQ*eq
            temp2 = prod_quat(Q,bar_quat(this.Qg)); % Q x Qg_bar
            temp3 = prod_quat(temp1,Q); % JlogQ*eq x Q 
            Q_dot = -1/2*prod_quat(temp2,temp3);
        end
        
        % Function that returns the time derivative of the orientation error.
        % Combination of eq. (20) and eq. (22)
        % @param[in] Q : Quaternion at current moment.
        % @param[in] Q_dot : Time derivative of the quaternion.
        function eq_dot = get_deq(this,Q,Q_dot)
            temp1 = prod_quat(this.Qg,bar_quat(Q));
            temp2 = prod_quat(Q_dot,bar_quat(Q));
            eq_dot = -2*quat_JQ(Q)*(prod_quat(temp1,prod_quat(temp1,temp2)));   
        end     
        
    end

    properties (Access = public)
        w % weights of dmp.
        clock % canonical clock of users choice.
        regressor % regressor object to obtain regressor values.
        az % Parameter for the dmp.
        bz % Parameter for the dmp.
        tau % Parameter for the dmp. 
        Q0 % Initial orientation.
        Q0d % Initial demonstrated orientation.
        Qg % Goal orientation.
        Qgd % Goal demonstrated orientation.
        S % Scaling for the dmp.
        Fx % Used instead of weigths.
    end

end


