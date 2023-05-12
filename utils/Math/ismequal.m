% Function to check equality for 2 values. Apparaently it is necessary.
function flag = ismequal(a,b,tol)
    if (nargin < 3)
        tol = 1e-6;
    end
    if length(a) > 1
        if length(a) == length(b)
            temp = a-b;
            n = length(a);
            flag = 1;
            for i=1:n
                if abs(temp(i)) > tol
                    flag = 0;
                end
            end
            return;
        elseif length(b) > 1
            error('Element 2 is a vector while Element 1 is an element');
        end
    end
    temp = a-b;
    if abs(temp) > tol
        flag = 0;
    else
        flag = 1;
    end
end