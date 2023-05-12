% Function that returns the time derivative of a quaternion. 
% The time derivative of the quaternion is calculated from the eq. (7) from Orientation_Exponential_Tracking.
% @param[in] Q : Quaternion from which the angle and the rotational axis will be extracted.
% @param[in] w : Angular velocity of the frame.
function  Q_dot = get_dotproduct(Q,w)
    J = get_JQ(Q);
    Q_dot = 1/2*J*w;
end