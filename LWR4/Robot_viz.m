clear all
close all


%% Create the KUKA LWR 4+ robot as a SerialLink object
robot = lwr_create();

%% Set the initial robot configuration
q = [0 pi/2 pi/3 -pi/2 0 0 0];

%% Move the robot base using the homogeneous transformation gw0. By default the robot base is g=I4
gw0 = [1 0 0 0.2;
       0 1 0 0.6;
       0 0 1 0.7;
       0 0 0 1];

robot.base = gw0;   
   
%% Plot the robot in the given configuration
robot.plot(q)