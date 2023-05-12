classdef dmp < handle
    methods (Access = public)
        % DMP constructor 
        %  @param[in] nBF: The number of basis functions(kernels).
        %  @param[in] T : Scale for canonical clock.
        %  @param[in] can_clock_index : Index to choose which canonical clock to use[0:exp,1:linear].
        function this = dmp(nBF,T,can_clock_index)
            if (nargin < 3)
                this.clock = canClock(T,1);
            else
                this.clock = canClock(T,can_clock_index);
            end
            this.w = zeros(nBF,1);
            this.regressor = regressor(nBF);
        end

        % Initialization function
        %  @param[in] T : Timesteps that correspond to y_data.
        %  @param[in] y_data : Demonstrated position trajectory.
        %  @param[in] y0 : Desired initial position
        %  @param[in] g : Desired final position(goal).
        function init(this,T,y_data,y0,g)
            this.trainDMP(T,y_data);
            this.y0d = y_data(:,1);
            this.gd = y_data(:,end);
            this.n_dof = length(y_data(:,1));
            if (nargin < 2)
                this.y0 = this.y0d;
                this.g = this.gd;
            else
                this.y0 = y0;
                this.g = g;
            end
            
            if this.n_dof == 3  
                % Get ks using : A novel...
                this.ksflag = 1;
            else
                % Get ks using : A reversible...
                this.ksflag = 0;
            end
            this.setks();
        end

        % Function that returns the reference trajectory.
        %  @param[in] t : Current time moment.
        function yx = getYx(this,t)
            x = this.clock.getPhase(t);
            phi = this.regressor.get_Fx(x);
            yx = this.ks*(this.w'*phi-this.y0d)+this.y0;
        end

        % Function that returns the reference velocity.
        %  @param[in] t : Current time moment.
        function yx_dot = getYx_dot(this,t)
            x = this.clock.getPhase(t);
            x_dot = this.clock.getPhaseDot(t);
            phi_dot = this.regressor.get_dFx(x,x_dot);
            yx_dot = this.ks*this.w'*phi_dot;
        end

        % Function that returns the reference acceleration.
        %  @param[in] t : Current time moment.
        function yx_ddot = getYx_ddot(this,t)
            x = this.clock.getPhase(t);
            x_dot = this.clock.getPhaseDot(t);
            x_ddot = this.clock.getPhaseDDot(t);
            phi_ddot = this.regressor.get_ddFx(x,x_dot,x_ddot);
            yx_ddot = this.ks*this.w'*phi_ddot;
        end


        % Function that simulates the system.
        %  @param[in] t : Current time moment.
        %  @param[in] y : Previous position.
        %  @param[in] y_dot : Previous velocity.
        %  @param[in] frep : Coupling term.
        %  @param[in] T : In case simulation time changes.
        function [y,dy,ddy] = simulation(this,t,y,y_dot,frep,T)
            if ( nargin < 5 )
                frep = 0;
            end
            if ( nargin == 6)
                this.set_clock(canClock(T));
            end
            % Global step of integration
            global dt
            % Desired stiffness parameter
            K = 400;
            % Desired damping parameter
            D = 40;

            % Get reference values
            yx = this.getYx(t);
            yx_dot = this.getYx_dot(t);
            yx_ddot = this.getYx_ddot(t);

            % Integration
            ddy = yx_ddot - D*(y_dot-yx_dot) - K*(y-yx)+frep;
            y = y + y_dot*dt;
            dy = y_dot+ddy*dt;
        end
        
        % Function that sets scaling method.
        %  @param[in] flag : Desired scale method (0 : A reversible, 1 : a
        %  novel).
        function set_scaleMethod(this,flag)
            if flag
                this.ksflag = 1;
                this.setks();
            else
                this.ksflag = 0;
                this.setks();
            end
        end

    end

    

    methods (Access = protected)
        % Scaling term setting function. The parameters below are parsed
        % with this pointer.
        %  If there are not 3 degrees of freedom the default is 1D scaling
        function setks(this)
            if  ( ~this.ksflag ) || ( this.n_dof ~= 3)
                this.ks = this.get1Dscaling();
            else
                this.ks = this.get3Dscaling();
            end
            % In order for ks to always be a matrix.
            [n,m] = size(this.getks());
            if n > m
                this.ks = diag(this.ks);
            end
        end
        
        % Setters.

        % Setter for weights. It is used in dmp_upd
        function set_w(this,W)
            this.w = W;
        end

        function set_g(this,g)
            this.g = g;
        end

        function set_y0(this,y0)
            this.y0 = y0;
        end
        
        function set_clock(this,clock)
            this.clock = clock;
        end

        % Getters.
        function g = get_g(this)
            g = this.g;
        end

        function gd = get_gd(this)
            gd = this.gd;
        end

        function y0 = get_y0(this)
            y0 = this.y0;
        end

        function y0d = get_y0d(this)
             y0d = this.y0d;
        end

        function regressor = get_regressor(this)
            regressor = this.regressor;
        end

        function clock = get_clock(this)
            clock = this.clock;
        end

        function n_dof = get_ndof(this)
            n_dof = this.n_dof;
        end

        function W = getW(this)
            W = this.w;
        end

        function ks = getks(this)
            ks = this.ks;
        end


    end

    methods (Access = private)
        % Train the DMP
        %  @param[in] T : Timesteps that correspond to y_data.
        %  @param[in] y_data : This DMP adapts based only on the position profile.
        function trainDMP(this,T,y_data)
            y_data = checkFormat(y_data);
            this.n_dof = y_data(:,1);
            this.y0d = y_data(:,1);
            C = this.clock.getPhase(T);
            psi = this.regressor.get_Fx(C);
            % Calculate Weights based on Least Squares TODO: IMPLEMENT LWR
            this.set_weights(psi,y_data);
        end
        
        % Function that sets the scaling term based on eq11 from this paper : A Reversible Dynamic Movement Primitive formulation
        % This function arguments are explained in scaling_term function.
        function ks = get1Dscaling(this)
            ks = (this.g-this.y0)./(this.gd-this.y0d);
        end

        % Function that sets the scaling term based on eq11 from this paper : A novel DMP formulation for global and frame independent spatial scaling in the task space.
        % This function arguments are explained in scaling_term function.
        % Function scaling should be in the same folder as dmp class
        function ks = get3Dscaling(this)
            [ks,~,~] = scaling(this.y0d,this.gd,this.y0,this.g);
        end

        % Function that adapts the weights of the dmp kernels based on least squares objective function.
        % @param[in] psi : Regressor values. This vector should be nBF x m where m is the total number of timestamps.
        % @param[in] yd : Demonstrated trajectory.
        % @param[out] w : Weights for the trajectory derived from Least Squares.
        % The weights should be : nBF x n_dof (for formula correspondance)
        function set_weights(this,psi,yd)
            this.w = ((yd*(psi'))/(psi*(psi')))';
        end
    end

    properties (Access = private)
        w % weights of dmp.
        ks % scaling term of reference.
        n_dof % degree of freedom.
        clock % canonical clock of users choice.
        regressor % regressor object to obtain regressor values.
        g % Desired goal.
        y0 % Desired initial position.
        y0d % Initial position extracted from demo.
        gd % Final position extracted from demo.
        ksflag % flag for ks.
    end

end


