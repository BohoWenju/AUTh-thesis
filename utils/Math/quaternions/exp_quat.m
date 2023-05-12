% Function that returns a quaternion based on a tangent vector x E R^3.
%  @param[in] x : Vector from which the quaternion will be retrieved.
function Q = exp_quat(x)
    n = norm(x);
    if ismequal(n,0) % tolerance = 1e-6
        Q = [1; 0; 0; 0];
        return;
    end
    Q = [cos(n); sin(n)*n*x];
end