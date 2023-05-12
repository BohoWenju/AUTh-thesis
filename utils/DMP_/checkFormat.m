% Change the format to correspond with the functions. 
% Should be nxm, where n is n_dof and m instances.
%   @param[in] x : input data.
function y = checkFormat(x)
    [n,m] = size(x);
    if n > m
        y = x';
        this.n_dof = m;
        return;
    end
    this.n_dof = n;
    y = x;
end