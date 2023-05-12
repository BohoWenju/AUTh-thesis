% Class to distinct the robots. Basically the controller will not be changed and when a new
% robot is to be included one has to only update the functions in this class by also adding the name of the robot
% in the definition.
classdef robot < handle
   methods (Access = public)

        % Constructor for class robot.
        %  @param[in] Hb : Homogeneous transform for the base of the robot.
        %  @param[in] robot_flag : Flag to indicate which robot is used.
        function this = robot(Hb, robot_flag)
            this.Hb = Hb;
            if ( nargin < 2)
                this.robot_flag = 1; % default is lwr4.
            else
                this.robot_flag = robot_flag;
            end
        end

        % Function to get the homogeneous transform of the robot's end-effector in respect to the base of the robot.
        % The particular usage of this function is to discriminate the robots here rather than in the main code.
        % Basically it solves the forward kinematics problem for a certain type of robot.
        %  @param[in] q : State vector that contains the angles of the robot's motors.
        function T = fkine(this, q)
            if this.robot_flag == 1
                T = lwr4_fkine(q);
            end
            T = this.Hb*T;
        end

        % Function to get the necessary matrices for the dynamics equation of the robot.
        %  @param[in] q : State vector that contains the angles of the robot's motors.
        %  @param[in] dq : State vector that contains the rotational speed/ velocity of the robot's motors.
        function [Je, M, C, G] = robotMatrices(this,q,dq)
            if this.robot_flag == 1
                if (nargin < 3)
                    Je = lwr4_jacob(q);
                    M = lwr4_inertia(q);
                    G = lwr4_gravity(q);
                    C = zeros(7,7);
                else
                    Je = lwr4_jacob(q);
                    M = lwr4_inertia(q);
                    G = lwr4_gravity(q);
                    C = lwr4_coriolis(q,dq);
                end
            end
        end

        % Function that returns the task space values based on the joint values.
        %  @param[in] q : State vector that contains the angles of the robot's motors.
        %  @param[in] dq : State vector that contains the rotational speed/ velocity of the robot's motors.       
        function [x,dx,Q,w] = joint2task(this,q,dq)
            T = this.fkine(q);
            [Je, ~, ~, ~] = this.robotMatrices(q);
            x = T(1:3,4);
            dx = Je*dq;
            Q = rotm2quat(T(1:3,1:3))';
            w = dx(4:6);
            dx = dx(1:3);
        end

        % Function that gets the orientation of the tool's end-effector based on the angles of the motors.
        % The output is in quaternion form and not in rotation matrix.
        % param[in] q : State vector that contains the angles of the robot's motors.
        function Q = get_quat(this,q)
            R = this.fkine(q);
            R = R(1:3,1:3);
            Q = rotm2quat(R);
            Q = make_unitary(Q);
        end

        function Hb = get_Hb(this)
            Hb = this.Hb;
        end
   end

   properties ( Access = private )
        Hb % Homogeneous transform of the robot's base.
        robot_flag % If 1 lwr4 if 0 anything else.
   end
end