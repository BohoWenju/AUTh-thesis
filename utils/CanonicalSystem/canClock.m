classdef canClock 
    methods (Access = public)
        % Canonical Clock constructor.
        % @param[in] scale : Parameter to define the scaling of the canonical system.
        % @param[in] can_clock_index : Flag to define whether the exponential canonical clock(0) should be used or the linear (~0)
        function this = canClock(scale,can_clock_index)
            if (nargin < 2)
                this.index = 1;
            else
                this.index = can_clock_index;
            end
            this.tau = scale;
            this.ac = 1;
        end
        
        %  @param[in] t : Moment on which the desired phase value should be calculated.
        function x = getPhase(this,t)
            if this.index
                x = (1/this.tau)*t;
            else
                x = exp(-(t*this.ac)/this.tau);
            end
        end

        %  @param[in] t : Moment on which the desired phase_dot value should be calculated.
        function x_dot = getPhaseDot(this,t)
            timesteps = length(t);
            if this.index
                x_dot = (1/this.tau)*ones(length(timesteps),1);
            else
                x_dot = -(this.ac/this.tau)*exp(-(t*this.ac)/(this.tau));
            end
        end

        %  @param[in] t : Moment on which the desired phase_ddot value should be calculated.
        function x_ddot = getPhaseDDot(this,t)
            timesteps = length(t);
            if this.index
                x_ddot = zeros(length(timesteps),1);
            else
                x_ddot = (this.ac^2/(this.tau^2))*exp(-(t*this.ac)/(this.tau));
            end
        end

        % Getter.
        function tau = get_tau(this)
            tau = this.tau;
        end

    end

    properties (Access = private)
        index % Flag to choose between desired canonical clock
        tau % Scaling parameter
        ac % Parameter for exponential canonical clock
    end
end