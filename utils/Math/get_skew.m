% Function that returns the skew matrix of a vector e.
%  @param[in] e : Vector e.
function a_skew = get_skew(a)
    a_skew = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end
