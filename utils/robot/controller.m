classdef controller < handle
    methods (Access = public)
        % Robot constructor 
        % Robot follows : M*d2q + C*dq + G = t + t_ext. t : is the torque the motor's produce.
        % t_ext : is the external torque applied to the robot at each timestep.
        % For the robot constructor:
        %  @param[in] robot : Class of robot. Until a more straightforward way is used a dictionary usage will
        %   be applied. 0: Whatever, 1: LWR4.
        %  The parameter size below depends on the type of control to be used. If joint space control is used then 
        %   the size must me n_dof x n_dof. Unless this is the case the size of position matrices must be 3 x 3 on task space. 
        %   The rotation matrices hold a single scalar value.
        %  @param[in] Kpos : Stiffness matrix for position control. 
        %  @param[in] Dpos : Damping matrix for position control.
        %  @param[in] Krot : Stiffness matrix for rotation control. This is not mandatory. 
        %  @param[in] Drot : Damping matrix for rotation control. This is not mandatory.
        function this = controller(robot,Kpos,Dpos,Krot,Drot)
            this.robot = robot;

            if ( nargin < 2 ) 
                this.Kpos = 400;
                this.Dpos = 2*sqrt(this.Kpos); % in order to be critically dumped.
                this.Krot = 400;
                this.Drot = 40;
                return;
            end
            this.Kpos = Kpos;
            this.Dpos = Dpos;
            this.Krot = Krot;
            this.Drot = Drot;
        end
        
        % Function that simulates the system.
        %  @param[in] t : Current time moment. Not needed but showcases simulation.
        %  @params[in] yd, yd_dot, yd_ddot : Position dmp output.
        %  @params[in] Qd, wd, wd_dot : Rotation dmp output. 
        %  @param[in] q : Current joing angles.
        %  @param[in] dq : Current angular velocity of joints.
        %  @param[in] fext : External GENERALIZED forces.
        function [q, dq, d2q, xe, Q] = simulation(this,t,xe,yd,yd_dot,yd_ddot,Qd,wd,wd_dot,q,dq,fext)
            if (nargin < 13)
                fext = zeros(6,1);
            end
            global dt 
            tau = 1;
            
            % Get matrices of robot.
            [Je, M, ~, ~] = this.robot.robotMatrices(q);

            % Get task space values.
            [x,dx,Q,w] = this.robot.joint2task(q,dq);
            
            % Inverse dynamics controller.
            up = this.pos_control(xe,dx,yd,yd_dot,yd_ddot);
%             [ur,~,~] = this.rot_control(Q,Qd,w,wd,wd_dot);
            temp = this.robot.fkine(q);
            x_axis = temp(1:3,1);
            ndy = make_unitary(yd_dot);
            theta = acos(x_axis'*ndy);
            ne = get_skew(x_axis)*ndy;
            eo = theta*ne;
            ur = -this.Drot*w-this.Krot*eo;
            u = [up; ur];

            % Closed loop inverse dynamics.
            c_torque = this.closed_loop(u,q,dq);
            ext_torque = (Je')*fext;

            % Integration.
            d2q = M\(c_torque)+M\ext_torque;
            dq = dq+d2q*tau*dt;
            q = q + dq*tau*dt;
            temp =lwr4_fkine(q);
            temp = this.robot.get_Hb()*temp;
            xe = temp(1:3,4);
        end

    end    
    methods (Access = private)

        %% Controller

        % Function that returns the necessary cartesian input to be parsed in the controller.
        %  @param[in] x : Position in space, x E R^3.
        %  @param[in] dx : Velocity in space, dx E R^3.
        %  @param[in] xd : Desired position in space.
        %  @param[in] xd_dot : Desired velocity in space.
        %  @param[in] xd_ddot : Desired acceleration in space.
        function [u, error] = pos_control(this,x,dx,xd,xd_dot,xd_ddot)
            % Get tracking error_dot
            error_dot = dx - xd_dot;
            % Get tracking error
            error = x - xd;
            u = xd_ddot-this.Kpos*error-this.Dpos*error_dot;
        end

        % Function that calculates the necessary input to be parsed in the controller of the robot
        % It takes into account a desired orientation trajectory and attempts to make the robot's end-effector
        % to follow this orientation. On this particular work as, as it is shown in the examples, the desired
        % orientation is calculated on each time-step and the desired angular velocity and angular acceleration
        % are set to 0. It is derived from Orientation_Exponential_Tracking. The particular methods are chosen because it
        % has been proven with experiments that they have the best results (NOTE DO THAT).
        % Quaternions are expected to be 4 x 1 and angular velocity/acceleration to be 3 x 1.
        % @param[in] Q : Current quaternion. The particular controller works in a feedback loop.
        % @param[in] Qd : Desired quaternion.
        % @param[in] w : Current angular velocity.
        % @param[in] wd : Desired angular velocity.
        % @param[in] wd_dot : Desired angular acceleration.
        % @param[out] u : Necessary torque for the robot's end effector.
        % @param[out] Qerror : Error quaternion. It is returned for validating.
        % @param[out] e : Norm of the logarithmic error.
        function [u,Qerror,e] = rot_control(this,Q,Qd,w,wd,wd_dot)

            if (nargin < 4)
                wd_dot = zeros(3,1);
                wd = zeros(3,1);
            end
            % use method i=l from leonidas koutras            
            % Getting quaternion error.
            Qerror = get_Qerror(Q,Qd);

            %     i=l eq. (21), (22)
            [theta,n] = quat2axis(Qerror);
            e = 2*log_quat(Qerror);
            Psie = theta*theta*1/2;
            % %     j=g
            R = rodriguez(n,theta);
            we = w-R*wd;
            dR = get_skew(we)*R*wd+R*wd_dot;
            u = dR-this.Drot*we-this.Krot*e;
            e = norm(e);
        end

        % Function that transforms the force and torque applied to the robot in task space
        %   to forces applied in the joint space.
        %  @param[in] u : Input to be parsed in the controller. This input must be mapped to it's joint space corresponding value.
        %  @param[in] q : State vector that contains the angles of the robot's motors.
        %  @param[in] dq : State vector that contains the rotational speed/ velocity of the robot's motors.
        function c_torque = closed_loop(this,u,q,dq)
            global dt % Get numerical step for integration.
            [Je, M, ~, ~] = this.robot.robotMatrices(q);
            qn = q + dq*dt;
            [Je_n,~,~,~] = this.robot.robotMatrices(q); 
%             Je_dot = (Je_n -  Je)/dt; % numerical differentiation. If a dJacob exists then substitute.
            Je_dot = zeros(6,7);
            L = Je*inv(M)*Je';
            L = inv(L); % Task space Inertia matrix.
            c_torque = (Je')*L*(u-Je_dot*dq); % Coriolis and gravity compensation has been applied.
        end        
        
    end
    
    properties (Access = private)
        Kpos % Position Stiffness.
        Dpos % Position Damping.
        Krot % Rotation Stiffness.
        Drot % Rotation Damping.
        robot % Robot obj.
    end
end

