% Function that returns the Jlog(Q) jacobian of a quaternion Q. eq. (19)
% @param[in] Q : Quaternion.
function J = getJlog(Q)
    [theta,n] = quat2axis(Q);
    theta = theta/2;
    if ismequal(theta,0)
        J = [zeros(1,3); eye(3)];
    else
        temp1 = -sin(theta)*n';
        temp2 = (sin(theta)/theta)*(eye(3)-n*n') + cos(theta)*(n*n');
        J = [temp1; temp2];
    end
end