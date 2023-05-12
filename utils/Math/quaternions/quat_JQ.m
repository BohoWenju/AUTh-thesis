% Function that returns the Jq(Q) jacobian of a quaternion Q. eq. (20)
% Note how it differs from get_JQ
% @param[in] Q : Quaternion.
function J = quat_JQ(Q)
    [theta,n] = quat2axis(Q);
    theta = theta/2;
    if ismequal(theta,0)
        J = [zeros(3,1) eye(3)];
    else
        s_theta = sin(theta);
        c_theta = cos(theta);
        temp = ((-s_theta+theta*c_theta)/(s_theta^2))*n;
        temp1 = (theta/s_theta)*eye(3);
        J = [temp temp1];
    end
end