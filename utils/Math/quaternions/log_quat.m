% Function that returns the logarithm of a quaternion.
% @param[in] Q : Quaternion from which the angle and the rotational axis will be extracted.
function logQ = log_quat(Q)
    [theta,n] = quat2axis(Q);
    if ismequal(n,zeros(3,1))
        logQ = [0 0 0]';
        return;
    end 
    logQ = n*(theta/2);
end
