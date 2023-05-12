% Must be in the same path with dmp class due to check format function.
classdef via < handle
    methods (Access = public)
        % Via point(s) constructor 
        %  @param[in] y : Desired position(s) for via points.
        %  @param[in] tvia : Desired moments on which the system must reach the via points.
        function this = via(y,tvia)
            y = checkFormat(y);
            n_via = length(y(1,:));
            nt_via = length(tvia(1,:));
            if (n_via ~= nt_via)
                error('Number of via points does not comply with number of 
                corresponding time moments\n');
            end
            this.n_via = n_via;
            this.y = y;
            this.tvia = tvia;
        end

        % Setter for via points. Has the same structure as constructor obviously.
        function set_via(this,y,tvia)
            y = checkFormat(y);
            n_via = length(y(1,:));
            nt_via = length(tvia(1,:));
            if (n_via ~= nt_via)
                error('Number of via points does not comply with number of 
                corresponding time moments.\n');
            end
            this.y = y;
            this.tvia = tvia;
        end

        % Getter for via points.
        %  @param[in] index : Desired via point to be extracted. Not necessary.
        function [via,tvia] = get_via(this,index)
            if (this.n_via < index )
                error('There is not a via points with this index.\n');
            end
            if (nargin < 1)
                via = this.y;
                tvia = this.tvia;
                return;
            end
            via = this.y(:,index);
            tvia = this.tvia(index);
        end

    end
    properties (Access = private)
        y % Position(s) of via point(s).
        tvia % Moment(s) on which the system must reach the via point(s).
        n_via % # of via points.
    end
end

