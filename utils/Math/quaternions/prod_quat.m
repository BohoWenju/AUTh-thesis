% Function that returns the product between 2 quaternions.
% @param[in] Qa : Quaternion A.
% @param[in] Qb : Quaternion B.
function Q = prod_quat(Qa,Qb)
    temp1 = bar_quat(Qa);
    temp2 = Qa(2:4);
    temp3 = Qa(1)*eye(3,3)+get_skew(temp2);
    temp = [temp1'; temp2 temp3];
    Q = temp*Qb;
%     Q = make_unitary(Q);
end