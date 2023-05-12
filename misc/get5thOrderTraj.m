% Function that generates 5th order trajectory with the constraints below:
% 1) Iniital position y0.
% 2) Final position yf.
% 3) Initial velocity 0.
% 4) Final velocity 0.
% 5) Initial acceleration 0.
% 6) Final acceleration 0.
% @param[in] y0 : Initial position.
% @param[in] yf : Final position.
% @param[in] Time : Complete set of time moments on which the trajectory is generated.
function [y, y_dot, y_ddot] = get5thOrderTraj(y0, yf, Time)

    y0 = y0(:);
    yf = yf(:);
    Time = Time(:)';
    
    T = Time(end);
    t = Time/T;
    
    y = y0 + (yf - y0) * (10*t.^3 - 15*t.^4 + 6*t.^5 );

    if (nargout > 1)
        y_dot = (yf - y0) * (30*t.^2 - 60*t.^3 + 30*t.^4 ) / T;
    end

    if (nargout > 2)
        y_ddot = (yf - y0) * (60*t - 180*t.^2 + 120*t.^3 ) / T^2;
    end

end