clear all
close all
setup();
%% =============== Fig7 Robot (Robot's end effector test) ===============
% This example simply uses the default dataset on each coordinate.
% The way this example is structured is : 
% 1. The forearm dataset is used to learn the dmp weights on a staight line in 3D space.
% 2. Place a robot in a certain poisition in space and move robot's end effector to hand and get it's orientation.
% 3. Initialize values for the simulation.
% 4. Simulate external forces applied to the robot by getting the forces in the form : f = k*z(z : closest point of arm - end_effector position)
% Calculate torque (torque = F x dr )and apply both the force and the torque to the 2 dmps as coupling terms in acceleration.
% 5. Plot the output in 3 positions.

%% =============== 1. ===============
%% Load Dataset
S = load('Forearm.mat');
% Extract Data
h = 0.2; % Height of the human arm.
for i = 1:3
    yd(i,:) = S.y(i,:);
    if i == 3
        yd(i,:) = yd(i,:) + h;
    end
    yd_dot(i,:) = S.y_dot(i,:);
    yd_ddot(i,:) = S.y_ddot(i,:);
end
Time = S.Time;
T = Time(end);
t0 = Time(1);
y0 = yd(:,1);
g = yd(:,end); 
gn = g;
global dt
dt = 2e-3;

%% DMP.
nBF = 30;
% Train the DMP
can_clock_index = 1;
% Choose optimization
optFlag = 0; % optimization based on acceleration.
dmp_up = dmp_upd(nBF,T,can_clock_index,optFlag);
dmp_up.init_upd(Time,yd,y0,g);
dmp_up.set_scaleMethod(1);
dmp_up.set_adapt(1);

% Create Constraint
Y = [y0 zeros(3,1) zeros(3,1)];
con = constraint(g,Y);

% Orientation dmp
rot_dmp = dmp_rot(nBF,T);

%% =============== 2. ===============
% Place the robot in the middle of the human arm and a little down.
R = eye(3);
base_pos = (yd(:,end) - yd(:,1))/2; % Place the robot in the middle.
base_pos(2) = base_pos(2) - 0.5;
base_pos(1) = base_pos(1) + 0.1;
Hb = [R base_pos; zeros(1,3) 1]; % Homogeneuous transform for the robot's base.
robot1 = lwr_create();
% Place the robot's end-effector to the human hand.
q0 = [0 pi/2 0 -pi/2 0 0 0]';
q1 = get_q(yd(:,1),q0,Hb);
q0 = q1;
% Construct robot and controller objects.
robot_obj = robot(Hb);
controller_obj = controller(robot_obj);
robot1.base = Hb;
% Get orientation of robot and use it to train the orientation dmp.
temp1 = robot_obj.fkine(q0);
Qd = rotm2quat(temp1(1:3,1:3))';
Qd = Qd.*ones(4,length(0:dt:T));
[eqd] = rot_dmp.init(Qd);

%% =============== 3. ===============


% Plan forearm movement. In this case it simply rotates around the hand.
phi = @(t)30*t/T; % Angle of rotation.

%% Set plots
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
title('$Tool-plot$','interpreter','latex', 'fontsize',13);
axis(ax, 'tight');
ax.XLim = ax.XLim + 0.05*(ax.XLim(2) - ax.XLim(1)) * [-1 1];
ax.YLim = [-0.05 1.4];
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off','location','southeast');
arm_pl = plot(nan, nan, 'Linestyle','-', 'linewidth',3, 'color',[0.85 0.33 0.1], 'DisplayName','Arm', 'Parent',ax);
dmp_pl = plot(nan, nan, 'Linestyle','-', 'linewidth',2, 'color',[0 0 1], 'DisplayName','Tool', 'Parent',ax);
yd_pl = plot(nan, nan, 'Linestyle',':', 'linewidth',2, 'color',[1 0 1],'DisplayName','DMP', 'Parent',ax);
uistack(dmp_pl, 'up', 10);
plot_count = 0;
plot_every = 10;
n_past = 7;

iters = 3;
timesteps = length(0:dt:T);
iter_plot_past = ceil([timesteps/4 timesteps/2 3*timesteps/4]);


%% =============== 4. ===============
i = 1;


% Make matrices.
timesteps = length(0:dt:T);
Q = zeros(4,timesteps);
w = zeros(3,timesteps);
w_dot = zeros(3,timesteps);
eq = zeros(3,timesteps);
eq_dot = zeros(3,timesteps);
eq_ddot = zeros(3,timesteps);
y = zeros(3,timesteps);
dy = zeros(3,timesteps);
ddy = zeros(3,timesteps);
q = zeros(7,timesteps);
dq = zeros(7,timesteps);
xe = zeros(3,timesteps);
% Initialize matrices.
Q(:,1) = Qd(:,1);
eq(:,i) = 2*log_quat(prod_quat(Qd(:,end),bar_quat(Qd(:,1))));
y(:,1) = temp1(1:3,4);
xe(:,1) = y(:,1);
Y = [y(:,i),dy(:,i),ddy(:,i)];
con.con_upd(g,Y);
q(:,1) = q0;

% Gains for the forces.
force_gain = 400;
torque_gain = 1;
% Simulation
for t = 0:dt:T
    % Rotate human arm.
    gn = rotz(phi(t))*g;
    temp = rotz(phi(t));
    Qd(:,i) = rotm2quat(temp);
    % Get force applied to the robot's end effector.
    [y1,~,f] = svc(y(:,i),y0,gn); % Get shortest vector to line.
    pos_coupling = (force_gain)*f;
    fmax = 0.5;
    a_ff = ( min([norm(pos_coupling)/ fmax, 1]) )^0.5;

    % Calculate torque.
    dr = make_unitary(y(:,i) - y(:,1)); % distance vector.
    rot_coupling = torque_gain*(get_skew(dr)*f);
    n_rot(i) = norm(rot_coupling);
    % Position DMP.
    [y(:,i+1),dy(:,i+1),ddy(:,i)] = dmp_up.simulation(...
    t,y(:,i),dy(:,i),1,con,pos_coupling,T,a_ff,0);
    
    % Rotation DMP.
    [Q(:,i+1), w(:,i+1),w_dot(:,i), eq(:,i+1), eq_dot(:,i+1), eq_ddot(:,i)] = rot_dmp.simulation(...
    t,eq(:,i),eq_dot(:,i),w(:,i),rot_coupling);

    % Put the DMP outputs in the controller.
    [q(:,i+1), dq(:,i+1), ddq(:,i), xe(:,i+1)] = controller_obj.simulation(...
    t,xe(:,i),y(:,i),dy(:,i),ddy(:,i),Q(:,i),w(:,i),w_dot(:,i),q(:,i),dq(:,i));
    
    % Translate the movement of the robot to the dmps.
    % Temp1 = robot_obj.fkine(q(:,i+1));
    % y(:,i+1) = Temp1(1:3,4);
    % [Je,~,~,~] = robot_obj.robotMatrices(q(:,i+1));
    % dx = Je*q(:,i);
    % w(:,i+1) = dx(4:6);
    % dy(:,i+1) = dx(1:3);
    
    
    Y = [y(:,i+1),dy(:,i+1),ddy(:,i)];
    con.con_upd(g,Y);
    i = i + 1;
    %% update plots
    dmp_pl.XData = [dmp_pl.XData xe(1,i-1)];
    dmp_pl.YData = [dmp_pl.YData xe(2,i-1)];
    
    arm_pl.XData = [y0(1) gn(1)];
    arm_pl.YData = [y0(2) gn(2)];
    
    if (~isempty(iter_plot_past) && i == iter_plot_past(1))
        past_arm_pl = copyobj(arm_pl, arm_pl.Parent);
        past_arm_pl.Color = [past_arm_pl.Color 0.5];
        set( past_arm_pl, 'HandleVisibility','off');
        iter_plot_past = iter_plot_past(2:end);
    end
    
    yd_pl.XData = [yd_pl.XData y(1,i)];
    yd_pl.YData = [yd_pl.YData y(2,i)];
    
    plot_count = plot_count + 1;
    if (plot_count == plot_every)
        plot_count = 0;
        drawnow();
        pause(0.001);
        %pause
    end
end


%% =============== 5. ===============
% Plots for orientation.
% In 3 distinctive time moments the orientation of the tool will be
% plotted along with the forearm.
Time = 0:dt:T;
time_moments = [ceil(timesteps/4) ceil(timesteps/2) ceil(3*timesteps/4)];
time_moments1 = [Time(time_moments(1)) Time(time_moments(2)) Time(time_moments(3))];
angle = [phi(time_moments1(1)) phi(time_moments1(2)) phi(time_moments1(3))];
for i = 1:3
    gn(:,i) = rotz(angle(i))*g;
end
s = Time/T;
k = y0 + s.*(gn(:,3)-y0);
figure();
title('3d-Arm and Robot plot');
hold on
plot3(k(1,:),k(2,:),k(3,:),'LineWidth',2, 'LineStyle','--','Color','blue','DisplayName','Arm Position')
hold on
robot1.plot(q(:,time_moments(3))');
xlabel("X");
ylabel("Y");
zlabel("Z");
hold off

% k = y0 + s.*(g-y0);
% figure();
% title('3d-Arm and Robot plot');
% hold on
% plot3(k(1,:),k(2,:),k(3,:),'LineWidth',2, 'LineStyle','--','Color','blue','DisplayName','Arm Position')
% hold on
% robot1.plot(q(:,1)');
% xlabel("X");
% ylabel("Y");
% zlabel("Z");
% hold off
for i = 1:timesteps
    eq(:,i) = 2*log_quat(get_Qerror(Qd(:,i),Q(:,i)));
end
Time = 0:dt:T+dt;
figure; hold on
subplot(3,1,1);
plot(Time,eq(1,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','Tool-epx');
grid on
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('eQx', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');
title("Quaternion erros");
subplot(3,1,2);
plot(Time,eq(2,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','Tool-epy');
grid on
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('eQy', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');
subplot(3,1,3);
plot(Time,eq(3,:),'LineWidth',2, 'LineStyle','-','Color','blue','DisplayName','Tool-epz');
grid on
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel('eQz', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');

Time = 0:dt:T;
legend_font = 17;
fig = figure;
ax_vec = [];
fig.Position(3:4) = [581 656];
% Plot Position
ax_vec = [ax_vec ax];
plot(Time,n_rot,'LineWidth',2, 'LineStyle','-','Color','red');
axis tight;
xlim([Time(1) Time(end)]);
xlabel('T [$s$]', 'interpreter','latex', 'fontsize',legend_font);
ylabel('Torque[$N*m$]', 'interpreter','latex', 'fontsize',legend_font);
hold off;

%% Functions.

function q = get_q(xe,q,g)
    dt = 2e-3;
    T = 4;
    Time = 0:dt:T;
    x = g*lwr4_fkine(q);
    x = x(1:3,4);
    [~, y_dot, ~] = get5thOrderTraj(x, xe, Time);
    i = 1;
    for t=0:dt:T
        J = lwr4_jacob(q);
        q_dot = pinv(J)*[y_dot(:,i); zeros(3,1)];
        q = q_dot*dt + q;
        i = i + 1;
    end
end

function setup()
    setPath();
end

function setPath()
    p = mfilename('fullpath');
    p = erase(p,'fig7_robot');
    addpath(genpath(p));
end
