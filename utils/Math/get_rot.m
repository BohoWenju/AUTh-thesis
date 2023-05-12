% Function that returns the rotational matrix who is derived from the vertical
% vector of 2 other vectors. Vectors must be unitary!
%  @param[in] a : Vector a.
%  @param[in] b : Vector b.
function k = get_rot(a,b)
    tol = 0.0001;
    if norm(a)  > 1+tol
        error('Vector a not unitary for rotation matrix');
    end
    if norm(b)  > 1+tol
        error('Vector b not unitary for rotation matrix');
    end
    temp = get_cross(a,b);
    if isnan(temp)
        if (a == -b)
            k = -eye(3);
            return;
        else
            k = eye(3);
            return;
        end
    end
    angle = acos((a')*b);
    k = rodriguez(temp,angle);
end



