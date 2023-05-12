% Function that extracts the angle and the rotational axis from a quaternion.
% @param[in] Q : Quaternion from which the angle and the rotational axis will be extracted.
% @param[in] tol : Slight tolerance. Not necessary.
function [theta,n] = quat2axis(Q,tol)
    if (nargin < 2)
        tol = 1e-8;
    end
    if ismequal(Q(1),1)
        theta = 2*atan2(0,1);
    else
        theta = 2*atan2(sqrt(1-Q(1)*Q(1)),Q(1));
    end
    if theta < tol
        theta = 0;
        n = [1 0 0]';
        return;
    elseif ismequal(theta,pi)
        theta = -pi;
    end
    n = sin(theta/2)*Q(2:4);
end