classdef constraint < handle
    methods (Access = public)
        % Constraint constructor
        %  @param[in] g : Desired target.
        %  @param[in] y0 : Desired initial position.
        %  @param[in] ypr : Previous state of the DMP.
        %  @param[in] via : Desired via point(s).
        function this = constraint(g,ypr,y0,via)
            if (nargin < 1)
                this.isGoal = 0;
                this.isY0 = 0;
                this.isYpr = 0;
                this.isVia = 0;
                return;
            end
            % y0, y_pr, via are on the example not implemented
            if (nargin < 2)
                set_goal(this,g);
                this.isY0 = 0;
                this.isYpr = 0;
                this.isVia = 0;
                return;
            end
            if (nargin < 3)
                set_goal(this,g);
                set_ypr(this,ypr);
                this.isY0 = 0;
                this.isVia = 0;
                return;
            end 
            if (nargin < 4)
                set_goal(this,g);
                set_y0(this,y0);
                set_ypr(this,ypr);
                this.isVia = 0;
                return;
            end 
            set_goal(this,g);
            set_y0(y0);
            set_ypr(ypr);
            set_via(via);
        end        

        % Function that updates the constraints based on new desired target,initial position
        % etc.
        % This function has the same input and structure as the constructor 
        % for obvious reasons.
        function con_upd(this,g,ypr,y0,via)
            % y0, y_pr, via are on the example not implemented
            if (nargin < 3)
                set_goal(this,g);
                this.isY0 = 0;
                this.isYpr = 0;
                this.isVia = 0;
                return;
            end
            if (nargin < 4)
                set_goal(this,g);
                set_ypr(this,ypr);
                this.isY0 = 0;
                this.isVia = 0;
                return;
            end 
            if (nargin < 5)
                set_goal(this,g);
                set_ypr(this,ypr);
                set_y0(this,y0);
                this.isVia = 0;
                return;
            end 
            set_goal(this,g);
            set_y0(this,y0);
            set_ypr(this,ypr);
            set_via(this,via);
        end

        % Deactivate function. Basically this function has the purpose to deactivate(hence the name)
        % certain constraints depending on user input. It is strutured that way for convenience.
        % If a constraint needs to remain active then the corresponding flag should be 1.
        %  @param[in] GoalFlag : Flag to show whether the goal constraint should remain.
        %  @param[in] Y0Flag : Flag to show whether the initial point constraint should remain.
        %  @param[in] YprFlag : Flag to show whether the previous state constraint should remain.
        %  @param[in] ViaFlag : Flag to show whether the via constraint should remain.
        function deactivate(this,GoalFlag,Y0Flag,YprFlag,ViaFlag)
            this.isGoal = GoalFlag;
            if ( nargin < 3 )
                return;
            else
                this.isY0 = Y0Flag;
                this.isYpr = YprFlag;
                this.isVia = ViaFlag;
            end
        end

        % Function to copy a constraint con.
        %  @param[in] con : Constraint to copy.
        function copy(this,con)
            if con.isGoal
                this.set_goal(con.get_goal);
            else
                this.isGoal = 0;
            end
            if con.isY0
                this.set_y0(con.get_y0);
            else
                this.isY0 = 0;
            end
            if con.isYpr
                this.set_ypr(con.get_ypr);
            else
                this.isYpr = 0;
            end
            if con.isVia 
                this.set_via(con.get_via);
            else
                this.isVia = 0;
            end  
        end

        % Function to remove Ypr from a constraint. This happens when this constraint
        % is to be removed.
        function remove(this)
            this.isYpr = 0;
        end
        
        % Getters.
        function g = get_goal(this)
            if this.isGoal
                g = this.g;
            else
                error('Goal is not parsed');
            end
        end
        
        function y0 = get_y0(this)
            if this.isY0
                y0 = this.y0;
            else
                error('Y0 is not parsed');
            end
        end
        
        function ypr = get_ypr(this)
            if this.isYpr
                ypr = this.ypr;
            else
                error('Ypr is not parsed');
            end
        end
        
        function via = get_via(this)
            if this.isVia
                via = this.via;
            else
                error('Via is not parsed');
            end
        end
        
        % Setters.
        %  @param[in] g : Desired target.
        function set_goal(this,g)
            this.g = g;
            this.isGoal = 1;
        end

        %  @param[in] y0 : Desired initial position.
        function set_y0(this,y0)
            this.y0 = y0;
            this.isY0= 1;
        end

        %  @param[in] ypr : Previous state of the DMP.
        function set_ypr(this,ypr)
            this.ypr = ypr;
            this.isYpr = 1;
        end

        %  @param[in] via : Desired via point(s).
        function set_via(this,via)
            this.via = via;
            this.isVia = 1;
        end

    end

    % Flags that are necessary.
    properties (Access = public)
        isGoal % Flag that states if goal constraint is active.
        isY0 % Flag that states if initial position constraint is active.
        isYpr % Flag that states if previous state constraint is active.
        isVia % Flag that states if via poitns constraint is active.
    end

    properties (Access = private)
        g % Desired goal.
        y0 % Desired initial position.
        ypr % Previous state of the DMP. Necessary for online updates.
        via % Desired via points.
    end
end


