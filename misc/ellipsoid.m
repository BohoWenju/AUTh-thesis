classdef ellipsoid < handle
    methods (Access = public)
        % Ellipsoid constructor 
        %  @param[in] c : Center of ellipsoid object.
        %  @param[in] S : Covariance of ellipsoid object.
        %  @param[in] d : Soft bounds of objects.
        function this = ellipsoid(c,S,k,d)
            if (nargin < 3)
                this.k = 0.1;
                this.d = 1;
            elseif (nargin < 4)
                this.k = k;
                this.d = 1;
            else
                this.k = k;
                this.d = d;
            end
            this.c = c;
            this.S = S*S;
            this.Si = inv(this.S);
            this.Sco = S;
        end

        % Function that returns the necessary force to be added in dmp as coupling term
        % so as to avoid hitting this object.
        function frep = repulsive_force(this,y)
            dV = this.calc_dV(y);
            frep = -this.k*dV;
        end

        % Getters.
        function c = get_c(this)
            c = this.c;
        end

        function d = get_d(this)
            d = this.d;
        end

        function S = get_S(this)
            S = this.S;
        end

        function Sco = get_Sco(this)
            Sco = this.Sco;
        end

    end


    methods (Access = private)

        % Get surface function of object.
        %  @param[in] y : Position in space. Depending on the output of this function
        % the y lies inside or outside of the object.
        function psi = surfun(this,y)
            psi = (((y-this.c)')*this.Si*(y-this.c)) - 1;
        end

        % Get derivative pf surface function of object.
        %  @param[in] y : Position in space. Depending on the output of this function
        % the y lies inside or outside of the object.
        function dpsi = dsurfun(this,y)
            dpsi = 2*this.Si*(y-this.c);
        end

        % Get error for object. Error is represented by a normalized distance of a position
        % y and the soft bounds d.
        %  @param[in] psi : Surface function result for a given position y.
        function e = surf_error(this,psi)
            e = ((psi - this.d)^2)/(this.d^2);
        end

        % Potential function derivative. Necessary to simulate repulsive forces when approaching object.
        % Potential function : V = -log(1-e).
        %  @param[in] y : Position in space.
        function dV = calc_dV(this,y)
            psi = this.surfun(y);
            if psi > this.d
                dV = zeros(length(y),1);
            else
                e = this.surf_error(psi);
                dpsi = this.dsurfun(y);
                temp = 2/(1-e);
                temp = temp*(sqrt(e)/this.d);
                dV = -temp*dpsi;
            end
        end

    end
    properties (Access = private)
        c % Center of object.
        S % Variance of object.
        Sco % Coovariance of object.
        Si % Inverse variance matrix of object.
        d % Soft bounds for object.
        k % Force gain.
    end
end

