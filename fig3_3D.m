clear all
close all
setup();
%% =============== Fig3 of DMP (3D) ===============
% This example simply uses the default dataset on each coordinate.
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
nBF = 25;
T = Time(end);
t0 = Time(1);
% One can decide to change goal or initial position i.e. g = yd(:,end) + 0.1
y0 = yd(:,1);
g = yd(:,end)+[0.2; 0.1; 0.1]; 

%% Train the DMP
can_clock_index = 1;
% Choose optimization
optFlag = 1;
dmp = dmp_upd(nBF,T,can_clock_index,optFlag);
dmp.init_upd(Time,yd,y0,g);
dmp.set_scaleMethod(1);

%% Create Constraint
Y = [y0 zeros(3,1) zeros(3,1)];
con = constraint(g,Y);



%% Simulate the DMP
% Set time step for integration
global dt
dt = 2e-3;

% Initial values for the integration
y = zeros(3,length(S.y));
dy = zeros(3,length(S.y));
ddy = zeros(3,length(S.y));
y(:,1) = y0;
dy(:,1) = zeros(3,1);
online = 1;
i=1;
for t=0:dt:T
    
    
    [y(:,i+1),dy(:,i+1),ddy(:,i)] = dmp.simulation(t,y(:,i),dy(:,i),online,con);
    % Get previous state.
    Y = [y(:,i+1),dy(:,i+1),ddy(:,i)];
    con.con_upd(g,Y);
    i = i + 1;
    
end
% New timesteps are required because of the initial values of the
% simulation.
simTime = 0:dt:T;
Timed = 0:dt:T+dt;


%% Plot results
ax_font = 13;
x_font = 16;
y_font = 16;
legend_font = 17;
% 3D plot
figure; hold on
plot3(y(1,:),y(2,:),y(3,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','DMP');
plot3(yd(1,:),yd(2,:),yd(3,:),'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo');
plot3(yd(1,1),yd(2,1),yd(3,1), 'LineStyle','None', 'Marker','o', 'Color',[0.5 1 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot3(yd(1,end),yd(2,end),yd(3,end),'LineStyle','None', 'Marker','x', 'Color',[1 0.5 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot3(y0(1),y0(2),y0(3), 'LineStyle','None', 'Marker','o', 'Color',[0 0.85 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot3(g(1),g(2),g(3),'LineStyle','None', 'Marker','x', 'Color',[0.85 0 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
grid on
axis equal
xlabel('X [$m$]', 'interpreter','latex', 'fontsize',legend_font);
ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',legend_font);
zlabel('Z [$m$]', 'interpreter','latex', 'fontsize',legend_font);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');





% Plot trajectories
for i=1:3
    fig = figure;
    ax_vec = [];
    fig.Position(3:4) = [581 656];
    % Plot Position
    ax = subplot(3,1,1); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot(Timed,y(i,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','DMP');
    plot(Time,yd(i,:),'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo');
    plot(t0,yd(i,1), 'LineStyle','None', 'Marker','o', 'Color',[0.5 1 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
    plot(T,yd(i,end),'LineStyle','None', 'Marker','x', 'Color',[1 0.5 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
    plot(t0,y0(i), 'LineStyle','None', 'Marker','o', 'Color',[0 0.85 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
    plot(T,g(i) ,'LineStyle','None', 'Marker','x', 'Color',[0.85 0 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off','DisplayName','$g$');
    ax.FontSize = ax_font;
    if i == 1
        ylabel('posx [$m$]', 'interpreter','latex', 'fontsize',y_font);
    elseif i == 2
        ylabel('posy [$m$]', 'interpreter','latex', 'fontsize',y_font);
    else
        ylabel('posz [$m$]', 'interpreter','latex', 'fontsize',y_font);
    end
    % title(['DoF $' num2str(i) '$'], 'interpreter','latex', 'fontsize',18);
    legend({}, 'interpreter','latex', 'fontsize',legend_font, 'Position',[0.2251 0.9396 0.5708 0.0425], 'Orientation','horizontal', 'Box','off');
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    % Plot Velocity
    ax = subplot(3,1,2); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Timed(1) Timed(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(Timed, dy(i,:), 'LineWidth',2.0, 'Color','blue');
    ax.FontSize = ax_font;
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',y_font);
    axis tight;
    xlim([Timed(1) Timed(end)]);
    hold off;
    
    % Plot Acceleration
    ax = subplot(3,1,3); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Timed(1) Timed(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(simTime, ddy(i,:), 'LineWidth',2.0, 'Color','blue');
    ax.FontSize = ax_font;
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',y_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',x_font);
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    linkaxes(ax_vec, 'x');
    for i=1:length(ax_vec)
       ax = ax_vec(i);
       ax.Box = 'on';
       ax.YLim = ax.YLim + 0.07*(ax.YLim(2)-ax.YLim(1))*[-1 1];
    end
    ax = ax_vec(i);
    ax.XLim(2) = ax.XLim(2) + 0.05;
end




function setup()
    setPath();
end

function setPath()
    p = mfilename('fullpath');
    p = erase(p,'fig3_3D');
    addpath(genpath(p));
end