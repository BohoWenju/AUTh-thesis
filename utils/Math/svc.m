% Function to find the point of a line that has the shortest vector between one point.
% This function does not use projection matrix but instead is derived from
% the dot product. Although one can be derived from another this is the product of my analysis.
% The necessary equations to prove the code below can be found in eq. 
% Note. To define a line in space one needs two points: point A and point B (usually).
% Note. If the point is not inside the field defined by the line then it returns either point A or point B.
%  @param[in] y : The point in space R^3.
%  @param[in] a : Point of line.
%  @param[in] b : Point of line.
%  @param[out] x : Point on line that is closest to the desired point. 
%  @param[out] dist : Distance between desired point and line.
%  @param[out] k : Vector from y to desired point.
% (HERE)
function [x, dist, k] = svc(y, a, b)

    % Unitary vector of line.
    n = getUnit(a,b);

    % Calculate angle between desired points and line points.
    if check_dot(getUnit(y, a), -n)
        if check_dot(getUnit(y, b), n)
            % The point is within line bounds.
            temp = y - a;
            temp1 = temp'*n;
            dist = sqrt(temp'*temp-temp1^2);
            x = a + n*temp1;
            k = x - y;
        else
            dist = norm(y-b);
            x = b;
            k = x - y;
        end
    else
        dist = norm(y-a);
        x = a;
        k = x - y;
    end
end

function n = getUnit(a,b)
    if norm(b-a) < 1e-7
        n = zeros(3,1);
        return;
    end
    n = (b-a)/norm(b-a);
end

function temp = dot_product(a,b)
    temp = a'*b;
end

function flag = check_dot(a,b)
    temp = dot_product(a,b);
    if temp > 0
        if temp <= 1
            flag = 1;
        else
            flag = 0;
        end
    else
        flag = 0;
    end
end