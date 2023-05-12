clear all
close all
setup();
%% =============== Fig5 of DMP (Coupling Arm force) ===============
% This example simply uses the default dataset on each coordinate.
% The way this example is structured is : 
% 1. The default dataset is used to learn the dmp weights on a staight line in 3D space.
% 2. The straight line then rotates by a small angle simulating the rotation of the arm.
% 3. On each moment t the position of dmp is calculated by also adding a coupling term that
% describes the forces acting on the end effector by the human arm.
% 4. Plotting the forces and the 3D output of the example.

%% Load Dataset
S = load('Data5th_order.mat');
% Extract Data
for i = 1:3
    yd(i,:) = S.y;
    yd_dot(i,:) = S.y_dot;
    yd_ddot(i,:) = S.y_ddot;
end
Time = S.Time;


%% Parameters for the DMP
nBF = 30;
T = Time(end);
t0 = Time(1);
% One can decide to change goal or initial position i.e. g = yd(:,end) + 0.1
yd(3,:) = zeros(1,length(yd(3,:)));
y0 = yd(:,1);
g = yd(:,end); 
gn = g;
%% Train the DMP
can_clock_index = 1;
% Choose optimization
optFlag = 0;
dmp_up = dmp_upd(nBF,T,can_clock_index,optFlag);
dmp_up.init_upd(Time,yd,y0,g);
dmp_up.set_scaleMethod(1);
dmp_up.set_adapt(1);

%% Create Constraint
Y = [y0 zeros(3,1) zeros(3,1)];
con = constraint(g,Y);

%% Simulate online movement of the arm
T = 10;
phi = @(t)30*t/T;

%% Initialization
% Set time step for integration
global dt
dt = 0.005;
% Initial values for the integration
timesteps = length(0:dt:T);
y = zeros(3,length(S.y));
dy = zeros(3,length(S.y));
ddy = zeros(3,length(S.y));
y(:,1) = y0;
dy(:,1) = zeros(3,1);

%% Set Plots
% path plot
figure; 
ax = axes();
hold(ax,'on'); grid on; box on;
plot(yd(1, :), yd(2, :), 'linewidth',2, 'color',0.4*[1 1 1], 'DisplayName','Demo', 'Parent',ax);
plot(y0(1), y0(2), 'Linestyle','none', 'Marker','o', 'linewidth',4, 'MarkerSize',10, 'color',[0 0.6 0], 'HandleVisibility','off', 'Parent',ax);
plot(g(1), g(2), 'Linestyle','none', 'Marker','x', 'linewidth',2, 'MarkerSize',13, 'color',[1 0.7 0.7], 'HandleVisibility','off', 'Parent',ax);
plot(gn(1), gn(2), 'Linestyle','none', 'Marker','x', 'linewidth',2, 'MarkerSize',13, 'color',[0.9 0 0], 'HandleVisibility','off', 'Parent',ax);
legend(ax, {}, 'fontsize',15, 'orientation','horizontal', 'Position',[0.2172 0.9355 0.5875 0.0624], 'box','off');
xlabel('X [$m$]', 'interpreter','latex', 'fontsize',14, 'Parent',ax);
ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',14, 'Parent',ax);
title('$DMP-accel$','interpreter','latex', 'fontsize',13);
axis(ax, 'tight');
ax.XLim = ax.XLim + 0.05*(ax.XLim(2) - ax.XLim(1)) * [-1 1];
ax.YLim = [-0.05 1.4];
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off','location','southeast');
arm_pl = plot(nan, nan, 'Linestyle','-', 'linewidth',3, 'color',[0.85 0.33 0.1], 'DisplayName','Arm', 'Parent',ax);
dmp_pl = plot(nan, nan, 'Linestyle','-', 'linewidth',2, 'color',[0 0 1], 'DisplayName','DMP', 'Parent',ax);
yd_pl = plot(nan, nan, 'Linestyle',':', 'linewidth',2, 'color',[1 0 1],'DisplayName','Desired', 'Parent',ax);
uistack(dmp_pl, 'up', 10);
plot_count = 0;
plot_every = 10;
n_past = 7;

iters = 3;
iter_plot_past = ceil([timesteps/4 timesteps/2 3*timesteps/4]);

% forces-stiffness plot
figure;
ax_f = {};
ax_f{1} = subplot(2,1,1);
f_plot = plot(nan, nan, 'linewidth',2, 'Color','magenta');
ylabel('[$N$]', 'interpreter','latex', 'fontsize',14);
title('$||F_{ext}||$', 'interpreter','latex', 'fontsize',14);
axis tight;
ax_f{2} = subplot(2,1,2);
a_ff_pl = plot(nan, nan, 'linewidth',2, 'color','green');
ylabel('$1-a_ff$', 'interpreter','latex', 'fontsize',14);
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',14);
axis tight;
for i=1:2
    ax_f{i}.XLim = [0 T];
end
ax_f{2}.YLim = [0 1];

%% Loop
online = 1;
i=1;
dt_camera = 0.06;
camera_counter = 0;
k = 400; % GAIN FOR ATTRACTION FORCES.
gn = rotz(10)*g;
for t = 0:dt:T
    %% Rotate arm continuously to simulate online movement.
    angle = phi(t);
    gn = rotz(angle)*g;
    
%     con.set_goal(gn);
    %% Camera input.
%     if ( camera_counter - dt_camera > 0 ) 
%         camera_counter = 0;
%         con.set_goal(gn);
%         flag = 1; % flag to know that camera has new frame.
%     end
%     camera_counter = camera_counter + dt;
%     
    %% Simulating forces
    [y1,~,f] = svc(y(:,i),y0,gn); % Get shortest vector to line.
    coupling_term = k*f;
    normforce(i) = norm(coupling_term); % Storing result for plot.
   
    fmax = 0.5;
    a_ff = ( min([norm(coupling_term)/ fmax, 1]) )^0.5;

    %% Simulation
    [y(:,i+1),dy(:,i+1),ddy(:,i)] = dmp_up.simulation(t,y(:,i),dy(:,i),online,con,coupling_term,T,a_ff);
    %     Get previous state.
    Y = [y(:,i+1),dy(:,i+1),ddy(:,i)];
    con.con_upd(g,Y);
    
    i = i + 1;
    
    %% update plots
    dmp_pl.XData = [dmp_pl.XData y(1,i-1)];
    dmp_pl.YData = [dmp_pl.YData y(2,i-1)];
    
    arm_pl.XData = [y0(1) gn(1)];
    arm_pl.YData = [y0(2) gn(2)];
    
    if (~isempty(iter_plot_past) && i == iter_plot_past(1))
        past_arm_pl = copyobj(arm_pl, arm_pl.Parent);
        past_arm_pl.Color = [past_arm_pl.Color 0.5];
        set( past_arm_pl, 'HandleVisibility','off');
        iter_plot_past = iter_plot_past(2:end);
    end
    
    yd_pl.XData = [yd_pl.XData y1(1)];
    yd_pl.YData = [yd_pl.YData y1(2)];
    
    f_plot.XData = [f_plot.XData t];
    f_plot.YData = [f_plot.YData norm(coupling_term)];
    
    a_ff_pl.XData = [a_ff_pl.XData t];
    a_ff_pl.YData = [a_ff_pl.YData 1-a_ff];
    
    plot_count = plot_count + 1;
    if (plot_count == plot_every)
        plot_count = 0;
        drawnow();
        pause(0.001);
        %pause
    end
    
end

%% Bin
% % New timesteps are required because of the initial values of the
% % simulation.
% simTime = 0:dt:T;
% Timed = 0:dt:T+dt;
% 
% %% Generate new line that occurs after rotating goal for the plots.
% yn = generate_line(y0,gn,Time);
% 
% %% Plot results. There is no need for 3D plot because z coordinates remain the same.
% ax_font = 13;
% x_font = 16;
% y_font = 16;
% legend_font = 17;
% % 3D plot
% figure; hold on
% plot(y(1,:),y(2,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','DMP');
% plot(yd(1,:),yd(2,:),'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo');
% plot(yn(1,:),yn(2,:),'LineWidth',2, 'LineStyle','-','Color','black','DisplayName','Desired');
% plot(yd(1,1),yd(2,1), 'LineStyle','None', 'Marker','o', 'Color',[0.5 1 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
% plot(yd(1,end),yd(2,end),'LineStyle','None', 'Marker','x', 'Color',[1 0.5 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
% plot(y0(1),y0(2), 'LineStyle','None', 'Marker','o', 'Color',[0 0.85 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
% plot(g(1),g(2),'LineStyle','None', 'Marker','x', 'Color',[0.85 0 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
% grid on
% axis equal
% xlabel('X [$m$]', 'interpreter','latex', 'fontsize',legend_font);
% ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',legend_font);
% legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off','location','northwest','Orientation','vertical');
% title('$DMP-accel$','interpreter','latex', 'fontsize',legend_font - 3);
% 


%% Functions
function y = generate_line(y0,g,Time)
    y = y0+(g-y0)*(Time/Time(end));
end


function setup()
    setPath();
end

function setPath()
    p = mfilename('fullpath');
    p = erase(p,'fig5_humancoupling');
    addpath(genpath(p));
end


function yd = get_nearest_pos_online(t, Tf, y0, g, y)
    % yd = argmin_{p} ||y - p||^2
    %      s.t. p = y0 + lambda*(g - y0)
    %           0 <= lambda <= 1
    g_y0 = g - y0;
    lambda = dot(g_y0, y-y0) / dot(g_y0, g_y0);
    if (lambda < 0), yd = y0;
    elseif (lambda > 1), yd = g;
    else, yd = y0 + lambda * g_y0;
    end
    
end