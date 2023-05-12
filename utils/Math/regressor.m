classdef regressor < handle
    methods (Access = public)
        % Regressor constructor 
        %  @param[in] nBF : Number of kernels.
        %  @param[in] c: The centers of the kernels.
        %  @param[in] h : Variance of basis functions.
        function this = regressor(nBF)
            this.nBF = nBF;
            this.set_centers();
            this.set_variance(1.5);
        end

        % Get Fx values
        %  @param[in] x : Canonical Clock output for a moment t.
        function Fx = get_Fx(this,x)
            if (x < 0)
                psi = zeros(length(this.nBF),1);
                psi(1) = 1;
            elseif (x >1)
                psi = zeros(length(this.nBF),1);
                psi(end) = 1;
            else
                psi = this.get_psi(x);
            end
            spsi = sum(psi,1);
            Fx = psi./spsi;
        end

        % Get dFx/dt values
        %  @param[in] x : Canonical Clock output for a moment t.
        %  @param[in] x_dot : Time derivative of Canonical Clock output for a moment t. 
        function dFx = get_dFx(this,x,x_dot)
            if (x < 0 || x > 1)
                dFx = zeros(length(this.nBF),1);
                return;
            end
            psi = this.get_psi(x);
            psi_dot = this.get_psiDot(x,x_dot);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);
        
            phi = psi / sum_psi; 
            dFx =  ( psi_dot - phi*sum_psi_dot ) / sum_psi; 
        end

        % Get ddFx/ddt values
        %  @param[in] x : Canonical Clock output for a moment t. 
        %  @param[in] x_dot : Time derivative of Canonical Clock output for a moment t. 
        %  @param[in] x_dot : Second Time derivative of Canonical Clock output for a moment t.
        function ddFx = get_ddFx(this,x,x_dot,x_ddot)
            if (x < 0 || x > 1)
                ddFx = zeros(length(this.nBF),1);
                return;
            end
            psi = this.get_psi(x);
            psi_dot = this.get_psiDot(x,x_dot);
            psi_ddot = this.get_psiDDot(x,x_dot,x_ddot);
            sum_psi = sum(psi);
            sum_psi_dot = sum(psi_dot);
            sum_psi_ddot = sum(psi_ddot);
        
            phi = psi / sum_psi; 
            phi_dot = ( psi_dot - phi*sum_psi_dot ) / sum_psi; 
            ddFx = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / sum_psi;   
        
        end

        % Getters.
        function c = get_c(this)
            c = this.c;
        end

        function h = get_h(this)
            h = this.h;
        end

        function nBF = get_nBF(this)
            nBF = this.nBF;
        end

    end

    methods ( Access = protected )
        % Get Psi values. Numerator of Fx components
        %  @param[in] x : Canonical Clock output for a moment t. 
        % The rest parameters are already explained.   
        function psi = get_psi(this,x)
            psi= exp(-this.h.*((x-this.c).^2));
        end

        % Get dPsi/dt values. d(Numerator of Fx components)/dt
        %  @param[in] x : Canonical Clock output for a moment t. 
        %  @param[in] x_dot : Time derivative of Canonical Clock output for a moment t. 
        % The rest parameters are already explained. 
        function psi_dot = get_psiDot(this,x,x_dot)
            psi = this.get_psi(x);
            a = (x-this.c)*x_dot;
            psi_dot = -2*this.h.*(psi.*a);
        end

        % Get ddPsi/ddt values. d(Numerator of Fx components)/dt
        %  @param[in] x : Canonical Clock output for a moment t. 
        %  @param[in] x_dot : Time derivative of Canonical Clock output for a moment t. 
        %  @param[in] x_dot : Second Time derivative of Canonical Clock output for a moment t.
        % The rest parameters are already explained. 
        function psi_ddot = get_psiDDot(this,x,x_dot,x_ddot)
            psi = this.get_psi(x);
            psi_dot = this.get_psiDot(x,x_dot);
            a = (x-this.c)*x_dot;
            a_dot = (x-this.c)*x_ddot + x_dot^2;
            psi_ddot= -2*this.h.*( psi_dot.*a + psi.*a_dot ); 
        end

        
    end

    methods ( Access = private )
        % Function that sets the centers of the kernels. Centers should be equally spaced in time between 0,1
        %   @param[in] nBF : # Basis Functions.
        function set_centers(this)
            this.c = linspace(0,1,this.nBF)';
        end

        % Function that sets the variance of the kernels. Variance should be determined on the previously set centers.
        % @param[in] ah : Scale parameter to determine the width of each basis function.
        function set_variance(this,ah)
            temp = this.c(2)-this.c(1);
            temp = (ah*temp)^2;
            temp = 1/temp;
            this.h = temp*ones(this.nBF,1);
        end
        
    end

    properties (Access = private)
        nBF % Number of kernels.
        c % The centers of the kernels.
        h % Variance of basis functions.
    end
end