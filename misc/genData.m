clear all
close all
setup();
%% =============== Generate 5th Order Trajectory ===============
% Creation of a dataset in order to use it as a demonstrated trajectory
% to use for the learning of a DMP.
% The particular dataset is derived from a 5th order trajectory.

% Integration step
dt = 1e-3; % 2e-2;

% Total duration
T = 8; % 3;

% Initial moment
t0 = 0;

%% Generate Trajectory

% Initial position
y0 = [0; 0; 0]; 

% Final position
yf = [0.25; .2; 0]; % 1;


Time = t0:dt:T;
[y, y_dot, y_ddot] = get5thOrderTraj(y0, yf, Time);

% Create the file with the dataset
filename = 'forearm.mat';
save(filename,'y','y_dot','y_ddot','Time');


%% Functions.
function setup()
    setPath();
end

function setPath()
    p = mfilename('fullpath');
    p = erase(p,'genData');
    addpath(genpath(p));
end

