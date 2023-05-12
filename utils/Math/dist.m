% Funciton that returns the distance between 2 vectors a, b.
%  @param[in] a : Vector a.
%  @param[in] b : Vector b.
function d = dist(a,b)
    d = sqrt(sum((a-b).^2));
end
