clear all
close all
setup();
%% =============== Fig2 of DMP (Online Upd (Both optimizations)) ===============

%% Load Dataset
S = load('Data5th_order.mat');
% Extract Data
yd = S.y;
yd_dot = S.y_dot;
yd_ddot = S.y_ddot;
Time = S.Time;

%% Parameters for the DMP
nBF = 25;
T = Time(end);
t0 = Time(1);
% One can decide to change goal or initial position i.e. g = yd(:,end) + 0.1
y0 = yd(:,1);
g = yd(:,end)+0.1; 

%% Train the DMP
can_clock_index = 1;
% Choose optimization
optFlag = 0;
dmp = dmp_upd(nBF,T,can_clock_index,optFlag);
dmp.init_upd(Time,yd,y0,g);


%% Create Constraint
Y = [y0 0 0];
con = constraint(g,Y);



%% Simulate the DMP
% Set time step for integration
global dt
dt = Time(2)-Time(1);

% Initial values for the integration
y(1) = y0;
dy(1) = 0;
online = 1;
i=1;
for t=0:dt:T
    
    
    [y(i+1),dy(i+1),ddy(i)] = dmp.simulation(t,y(i),dy(i),1,con);
    con.con_upd(g,[y(i+1),dy(i+1),ddy(i)]);
    i = i + 1;
    
end
% New timesteps are required because of the initial values of the
% simulation.
Timed = 0:dt:T+dt;
simTime = 0:dt:T;

%% Plot results
figure; hold on
plot(Timed,y,'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','DMP');
plot(Time,yd,'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo');
plot(t0,yd(:,1), 'LineStyle','None', 'Marker','o', 'Color',[0.5 1 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot(T,yd(:,end),'LineStyle','None', 'Marker','x', 'Color',[1 0.5 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot(t0,y0, 'LineStyle','None', 'Marker','o', 'Color',[0 0.85 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot(T,g ,'LineStyle','None', 'Marker','x', 'Color',[0.85 0 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
grid on
axis equal
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');

function setup()
    setPath();
end

function setPath()
    p = mfilename('fullpath');
    p = erase(p,'fig2_updmp');
    addpath(genpath(p));
end