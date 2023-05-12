% Function to get unit vector between two points a, b.
%  @param[in] a : Point a.
%  @param[in] b : Point b.
function u = unit_vector(a,b)
    u = a-b;
    u = u./dist(a,b);
end
