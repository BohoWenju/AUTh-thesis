clear all
close all
setup();
%% =============== Fig6 Orientation ( Orientation DMP test) ===============
% This example simply uses the default dataset on each coordinate.
% The way this example is structured is : 
% 1. Load the dataset/
% 2. A rotation dmp will be trained using the quaternion dataset from step 2.
% 3. A new orientation Initial position and Goal will be chosen.
% 4. The dynamical system presented in the work : A correct formulation for the Orientation Movement Primitives for robot control in the Cartesian Space 
% will be simulated.
% 5. Finally, the plots for the orientation error will be demonstrated.

%% =============== 1. ===============
S = load('ori_data.mat');
Qd = S.Q_data;
T = S.t_data(end);
global dt
dt = S.t_data(2)-S.t_data(1);



%% =============== 2. ===============
% First create the rotation dmp obj.
nBF = 60;
dmp = dmp_rot(nBF,S.t_data(end));

% Initialize and train the DMP.
[eqd] = dmp.init(Qd,S.t_data);

%% =============== 3. ===============
% goal_angle = qd(:,end);
% Qg = rotm2quat(rotz(goal_angle))';
Qg = Qd(:,end);

% Assume for now that no change will occur in the initial orientation.

%% =============== 4. ===============
% Initial conditions.
i = 1;

% Make matrices.
timesteps = length(0:dt:T);
Q = zeros(4,timesteps);
w = zeros(3,timesteps);
w_dot = zeros(3,timesteps);
eq = zeros(3,timesteps);
eq_dot = zeros(3,timesteps);
eq_ddot = zeros(3,timesteps);

% Initialize matrices.
eq(:,i) = 2*log_quat(prod_quat(Qg,bar_quat(Qd(:,1))));
eq_dot(:,i) = zeros(3,1);
Q(:,1) = Qd(:,1);
w(:,i) = zeros(3,1);

for t = 0:dt:T
    [Q(:,i+1), w(:,i+1), w_dot(:,i), eq(:,i+1), eq_dot(:,i+1), eq_ddot(:,i)] = dmp.simulation(t,eq(:,i),eq_dot(:,i),w(:,i));
    i = i + 1;
end

%% =============== 5. ===============
% Plot results
% One also needs the error during the demonstration.
Timed = 0:dt:T+dt;
Time = S.t_data;
figure; hold on
subplot(3,1,1);
plot(Timed,eq(1,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','DMP-epx');
hold on
plot(Time,eqd(1,:),'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo-epx');
grid on
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('eQx', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');
subplot(3,1,2);
plot(Timed,eq(2,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','DMP-epy');
hold on
plot(Time,eqd(2,:),'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo-epy');
grid on
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('eQy', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');
subplot(3,1,3);
plot(Timed,eq(3,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','DMP-epz');
hold on
plot(Time,eqd(3,:),'LineWidth',2, 'LineStyle','--','Color','green','DisplayName','Demo-epz');
grid on
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('eQz', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');






%% Functions.
function setup()
    setPath();
end

function setPath()
    p = mfilename('fullpath');
    p = erase(p,'fig6_orientation');
    addpath(genpath(p));
end
