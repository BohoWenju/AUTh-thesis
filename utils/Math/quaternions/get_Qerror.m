% Function that returns the error between two orientations in space 
% if they are described in quaternion form
%  @param[in] Qa : First quaternion.
%  @param[in] Qb : Second quaternion.
function Qe = get_Qerror(Qa,Qb)
    temp1 = bar_quat(Qa);
    temp2 = Qa(2:4);
    temp3 = Qa(1)*eye(3,3)+get_skew(temp2);
    temp = [temp1'; temp2 temp3];
    Qe = temp*bar_quat(Qb);
    Qe = make_unitary(Qe);
    % Constrain the scalar part of the quaternion in order to constrain the error
    % angles in the range [-pi, pi).
    Qe(1) = abs(Qe(1));
end