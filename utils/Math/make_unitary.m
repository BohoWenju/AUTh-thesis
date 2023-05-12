% Function that transforms a vector no matter the size to a unitary vector.
% @param[in] a : The vector to be transformed.
function u = make_unitary(a)
    temp = norm(a);
    if temp <= 1e-9
        u = a;
    else
        u = a ./ temp;
    end
end