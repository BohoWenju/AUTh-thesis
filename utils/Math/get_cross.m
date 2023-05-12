% Function that gets the unitary cross product of 2 vectors. It uses the skew matrix of one 
% and the equivalency : a x b = skew(a)*b
%  @param[in] a : Vector a.
%  @param[in] b : Vector b.
function temp = get_cross(a,b)
    temp = get_skew(a)*b;
    temp = make_unitary(temp);
end

