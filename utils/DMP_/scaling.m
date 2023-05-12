% Function that sets the scaling term based on eq11 from this paper : A novel DMP formulation for global and frame independent spatial scaling in the task space.
%  @param[in] y0d: Demonstrated initial position.
%  @param[in] gd : Demonstrated final position.
%  @param[in] y0 : Desired initial position.
%  @param[in] g : Desired final position.
function [ks,sg,Rg] = scaling(y0d,gd,y0,g)
    nd = unit_vector(gd,y0d);
    n = unit_vector(g,y0);

    % magnitude
    temp = dist(g,y0);
    sg = dist(gd,y0d);
    temp = temp/sg;
    % rotation
    Rg = get_rot(nd,n);
    ks = temp*Rg;

end
