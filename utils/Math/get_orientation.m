% This function returns the necessary orientation relative to the previous one.
% in order for the robot's end effector to always remain vertical to the arm.
% The first orientation is set by default as vertical so it can be set as an initial condition.
% The output of this function is sent to the rotation controller.
%  @param[in] t : Only needed to distribute the initialization from the other moments.
%  @param[in] dy : previous velocity. Translational velocity.
%  @param[in] ndy : current velocity. Translational velocity.
% It is assumed that the 2 vectors of velocity as well as the quaternions are vectors descirbed as columns.
% For now this function works only in 2D.
function Q = get_orientation(dy,ndy,t)

    persistent k0 dy0;
    
    % First make the velocity vectors unitary.
    dy = make_unitary(dy);
    if ismequal(dy,zeros(3,1))
        dy = dy0;
    end
    ndy = make_unitary(ndy);
    
    % In inititalization dy, ndy are the points of the hand and the elbow(or elbow and shoulder).
    if (nargin < 3)
        % Slight tolerarnce is included because atan2 gets only real inputs.
        tol = 1e-8;
        % Get the angle between the two vectors.
        cos_angle = (dy')*ndy;
        angle = atan2(sqrt(1+tol-cos_angle^2),cos_angle); % This angle is also the angle between the two vertical axes to the motion.

        k = get_cross(dy,ndy);

        % If the 2 vectors are parallel to each other then the rotation axis cannot be found.
        % So in that case the previous axis is held and used.
        if isnan(k)
            k = k0;
        else
            if ismequal(k,zeros(3,1))
                k = k0;
            else
                k0 = k;
            end
        end
        % The second orientation can be found by rotating the end effector by the angle above around the first rotation axis, which is derived from the previous quaternion.
        R = rodriguez(k,angle);
        Q = rotm2quat(R);
    else
%         k0 = get_cross(dy,ndy);
        k0 = [0;0;1];
        dy0 = ndy;
        Q = zeros(4,1);
        return;
    end
    

end

