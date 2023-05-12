% Function that returns the matrix J(Q) of a quaternion.
% This matrix is needed for the calculation of the time derivative of a quaternion.
%  @param[in] Q : Desired quaternion.
function JQ = get_JQ(Q)
    JQ = [ -Q(2:4)'; Q(1)*eye(3)-get_skew(Q(2:4))];
end