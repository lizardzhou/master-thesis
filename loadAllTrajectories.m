%% Script for loading data of all particle trajectories from importFluentData.m
% When the information of particle trajectories is too large to export to one file, 
% it should be splitted into 3 parts for faster processing. Time
% consumption for loading each part is recorded.
% The generated cell "traj" contains info of all single particle trajectories.
% Information of each trajectory is stored in each element of the cell.
% Original trajectory files should be loaded with ascending diameters.
% Original trajectory files should contain information in the way described below:
%   (:,1) x
%   (:,2) y
%   (:,3) z
%   (:,4) residence time
%   (:,5) particle diameter calculated by Fluent
%   (:,6) particle diameter calculated by UDF drying kinetics
%   (:,7) moisture content
%   (:,8) water fraction calculated by Fluent
%   (:,9) relative drying rate
%   (:,10) single particle mass
%   (:,11) particle temperature
%   (:,12) particle u_x
%   (:,13) particle u_y
%   (:,14) particle u_z
%   (:,15) gas u_x
%   (:,16) gas u_y
%   (:,17) gas u_z
%   (:,18) time step (cumulative sum of time step is residence time)

% by Xiye Zhou, Oct. 2021
%% Load data part 1
fprintf('Time for loading diameter 1 - 22 \x03bcm: \n'); % \x03bc: unicode of Greek letter mu
tic
rawData1 = readtable('simulationData/150ksteps-experimentBC-res0.99/d01-22loaded.txt');
data1 = table2array(rawData1);
clear rawData1; % save memory
toc
%% Load data part 2
fprintf('Time for loading diameter 27 - 42 \x03bcm: \n');
tic
rawData2 = readtable('simulationData/150ksteps-experimentBC-res0.99/d27-42loaded.txt');
data2 = table2array(rawData2);
clear rawData2; % save memory
toc
%% Load data part 3
fprintf('Time for loading diameter 47 - 60 \x03bcm: \n');
tic
rawData3 = readtable('simulationData/150ksteps-experimentBC-res0.99/d47-60loaded.txt');
data3 = table2array(rawData3);
clear rawData3; % save memory
toc
%% sort out single trajectories for all data parts
traj1 = sortSingleTraj(data1);
traj2 = sortSingleTraj(data2);
traj3 = sortSingleTraj(data3);
%% Merge all trajectories in cell traj
traj = [traj1;traj2;traj3];
fprintf('All single trajectories are stored in cell "traj"! \n');
%% Count # of incomplete trajectories
maxNumIter = 150e3;
c = 0;
for i = 1:size(traj,1)
    if size(traj{i,1},1)  > maxNumIter
        c = c + 1;
    end
end
fprintf('%g trajectories are tracked. According to Fluent: %g complete, %g incomplete. \n \n', size(traj,1), size(traj,1)-c, c);
clear maxNumIter c