% Returns a rotational matrix based on rodriguez formula for 2 vectors.
%  @param[in] k : Vector on which the rotation occurs.
%  @param[in] angle : Angle for rotation.
function R = rodriguez(k,angle)
    ske = get_skew(k);
    R = eye(3,3)+sin(angle)*ske+ske*ske*(1-cos(angle));
end