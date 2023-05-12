clear all
close all
setup();
%% =============== Fig4 of DMP (Coupling repulsive forces) ===============
% This example simply uses the default dataset on each coordinate.
% The way this example is structured is : 
% 1. The default dataset is used to learn the dmp weights on a staight line in 3D space.
% 2. Objects then are inserted in the line to block the dmp from following the path.
% 3. On each moment t the position of dmp is calculated by also adding a coupling term that
% describes the repulsive forces acting on the end effector.
% 4. Plotting the forces and the 3D output of the example.

%% Load Dataset
S = load('Data5th_orderCoupling.mat');
% Extract Data
yd = S.y;
yd_dot = S.y_dot;
yd_ddot= S.y_ddot;
Time = S.Time;

%% Create two ellipsoid obstacles
c1 = [0.45;0.1; 0];
S1 = diag([0.2; 0.2; 0.2]);
obj1 = ellipsoid(c1,S1,1,1);

c2 = [1.25; 0.8; 0];
S2 = diag([0.3; 0.15; 0.2]);
obj2 = ellipsoid(c2,S2,1,1);



%% Parameters for the DMP
nBF = 25;
% One can decide to change goal or initial position i.e. g = yd(:,end) + 0.1
y0 = yd(:,1);
g = yd(:,end); 
T = Time(end);
t0 = Time(1);


%% Train the DMP
can_clock_index = 1;
% Choose optimization
optFlag = 1;
dmp_up = dmp_upd(nBF,T,can_clock_index,optFlag);
dmp_up.init_upd(Time,yd,y0,g);
dmp_up.set_scaleMethod(1);

% Use typical dmp
dmp_og = dmp(nBF,T,can_clock_index);
dmp_og.init(Time,yd,y0,g);

%% Create Constraint
Y = [y0 zeros(3,1) zeros(3,1)];
con = constraint(g,Y);


%% Simulate the DMP
% Set time step for integration
global dt
dt = 1e-3;
% Initial values for the integration
% For updated dmp.
y = zeros(3,length(S.y));
dy = zeros(3,length(S.y));
ddy = zeros(3,length(S.y));
y(:,1) = y0;
dy(:,1) = zeros(3,1);

% For typical dmp.
y1 = zeros(3,length(S.y));
dy1 = zeros(3,length(S.y));
ddy1 = zeros(3,length(S.y));
y1(:,1) = y0;
dy1(:,1) = zeros(3,1);

online = 1;
i=1;
frep = zeros(3,1);
for t=0:dt:T
    
    % Updated DMP.
    % Get repulsive forces to avoid obstacles.
    frep = frep + obj1.repulsive_force(y(:,i));
    frep = frep + obj2.repulsive_force(y(:,i));
    coupling_term = frep;
    normforce(i) = norm(coupling_term);
    
    % Simulation.
    [y(:,i+1),dy(:,i+1),ddy(:,i)] = dmp_up.simulation(t,y(:,i),dy(:,i),online,con,coupling_term);
    % Dmp online update of constraints.
    Y = [y(:,i+1),dy(:,i+1),ddy(:,i)];
    con.con_upd(g,Y);
    
    % Typical DMP.
    % Get repulsive forces to avoid obstacles.
    frep = zeros(3,1);
    frep = frep + obj1.repulsive_force(y1(:,i));
    frep = frep + obj2.repulsive_force(y1(:,i));
    coupling_term = frep;
    % Simulation.
    [y1(:,i+1),dy1(:,i+1),ddy1(:,i)] = dmp_og.simulation(t,y1(:,i),dy1(:,i),coupling_term);
    normforce1(i) = norm(coupling_term);
    
    
    i = i + 1;
    frep = zeros(3,1);
    
end
% New timesteps are required because of the initial values of the
% simulation.
simTime = 0:dt:T;
Timed = 0:dt:T+dt;

%% Generate new line that occurs after rotating goal for the plots.
[h_obj1,s_obj1] = generate_object(obj1);
[h_obj2,s_obj2] = generate_object(obj2);

%% Plot results. There is no need for 3D plot because z coordinates remain the same.
ax_font = 13;
x_font = 16;
y_font = 16;
legend_font = 17;
% 3D plot
figure; hold on
plot(y(1,:),y(2,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','UpDMP');
plot(y1(1,:),y1(2,:),'LineWidth',2, 'LineStyle','--','Color','magenta','DisplayName','DMP');
plot(yd(1,:),yd(2,:),'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo');
plot(h_obj1(1,:),h_obj1(2,:),'LineWidth',2, 'LineStyle','-','Color','black','HandleVisibility','off');
plot(s_obj1(1,:),s_obj1(2,:),'LineWidth',2, 'LineStyle','-','Color',[0.9 0.6 0.4],'HandleVisibility','off');
plot(h_obj2(1,:),h_obj2(2,:),'LineWidth',2, 'LineStyle','-','Color','black','HandleVisibility','off');
plot(s_obj2(1,:),s_obj2(2,:),'LineWidth',2, 'LineStyle','-','Color',[0.9 0.6 0.4],'HandleVisibility','off');
plot(yd(1,1),yd(2,1), 'LineStyle','None', 'Marker','o', 'Color',[0.5 1 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot(yd(1,end),yd(2,end),'LineStyle','None', 'Marker','x', 'Color',[1 0.5 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot(y0(1),y0(2), 'LineStyle','None', 'Marker','o', 'Color',[0 0.85 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot(g(1),g(2),'LineStyle','None', 'Marker','x', 'Color',[0.85 0 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
grid on
axis equal
xlabel('X [$m$]', 'interpreter','latex', 'fontsize',legend_font);
ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',legend_font);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');


% Plot forces.
figure; hold on
plot(simTime,normforce,'LineWidth',2, 'LineStyle','-','Color','red','DisplayName','UpForce');
plot(Time,normforce1,'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','DMPForce');
grid on
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',legend_font);
ylabel('F [$N$]', 'interpreter','latex', 'fontsize',legend_font);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');


% Plot trajectories
for i=1:2
    fig = figure;
    ax_vec = [];
    fig.Position(3:4) = [581 656];
    % Plot Position
    ax = subplot(3,1,1); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot(Timed,y(i,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','UpDMP');
    plot(Timed,y1(i,:),'LineWidth',2, 'LineStyle','--','Color','magenta','DisplayName','DMP');
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
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(Timed, dy(i,:), 'LineWidth',2.0,'LineStyle','-', 'Color','blue','DisplayName','UpDMP');
    plot(Timed, dy1(i,:), 'LineWidth',2.0,'LineStyle','--', 'Color','magenta','DisplayName','DMP');
    ax.FontSize = ax_font;
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',y_font);
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    % Plot Acceleration
    ax = subplot(3,1,3); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(simTime, ddy(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(simTime, ddy1(i,:), 'LineWidth',2.0,'LineStyle','--', 'Color','magenta','DisplayName','DMP');
    ax.FontSize = ax_font;
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',y_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',x_font);
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    linkaxes(ax_vec, 'x');
    for j = 1:length(ax_vec)
       ax = ax_vec(j);
       ax.Box = 'on';
       ax.YLim = ax.YLim + 0.07*(ax.YLim(2)-ax.YLim(1))*[-1 1];
    end
    ax = ax_vec(i);
    ax.XLim(2) = ax.XLim(2) + 0.05;
end


function [y, y1] = generate_object(obj)

    
    S1 = obj.get_Sco();
    S = obj.get_S()+obj.get_S()*obj.get_d();
    S2 = sqrt(S);
    theta = linspace(0,2*pi,200);
    temp = S1*[cos(theta);sin(theta);zeros(1,200)];
    y = obj.get_c() + temp; % hard bounds
    temp = S2*[cos(theta);sin(theta);zeros(1,200)];
    y1 = obj.get_c() + temp; % soft bounds
    
end


function setup()
    setPath();
end

function setPath()
    p = mfilename('fullpath');
    p = erase(p,'fig4_coupling');
    addpath(genpath(p));
end