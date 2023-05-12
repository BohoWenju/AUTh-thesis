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
function [u,Qerror,e] = rot_control(Q,Qd,w,wd,wd_dot)

    if (nargin < 4)
        wd_dot = zeros(3,1);
        wd = zeros(3,1);
    end
    % use method i=l from leonidas koutras
    % this will be changed when the particular function is applied to the robot.
    % Damping and Stiffness parameters.
    d=8;
    k=16;
    
    % Getting quaternion error.
    Qerror = get_error(Q,Qd);

%     i=l eq. (21), (22)
    [theta,~] = rot_axis(Qerror);
    e = 2*log_quat(Qerror,theta);
    Psie = theta*theta*1/2;
% %     j=g
    R = rodriguez(Qerror);
    we = w-R*wd;
    dR = get_skew(we)*R*wd+R*wd_dot;
    u = dR-d*we-k*e;
    e = get_norm(e);
end


