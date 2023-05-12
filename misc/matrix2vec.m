% Get vector from matrix.
% This function returns all the elements of a matrix in a vector
%  @param[in] M : Matrix.
function v = matrix2vec(M)
    [~,m] = size(M);
    k = 1;
    v = zeros(m,1);
    for i = 1:m
        for j = 1:3
            if M(j,i) ~= 0
                v(k) = M(j,i);
                k = k + 1;
            end
        end
    end
end