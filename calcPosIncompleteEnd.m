function [distanceWall, distanceWallArray, numInTraj] =calcPosIncompleteEnd(d, maxNumIter, dInjection, traj_dInjection)
% Calculates the distance to chamber wall of the end point of incomplete
% trajectories for a given injection diameter d
% Input:
%   d: user-input injection diameter in meter
%   maxNumIter: maximal number of data points per trajectory
%   dInjection: lists all injection diameters ascendingly. Generated in
%   dataAvgHeight.m by sortDiamTraj.m
%   traj_dInjection: cell of size numClass*1. Each cell element contains
%   particle of same injection diameter. Generated in dataAvg.m by sortDiamTraj.m
% Return:
%   distanceWallArray: calculation of wall distance
%       (:,1) trajectory ID in Fluent
%       (:,2) x of end point
%       (:,3) y of end point
%       (:,4) z of end point
%       (:,5) radius of wall
%       (:,6) distance of end point to wall
%   trajNum: trajectory number in cell "traj". The trajectory can be treated as
%   complete one as coarse product because of the very low height of end point 

% by Xiye Zhou, Oct. 2021
%% Geometry of drying chamber
geo_rUpper = 0.397; % radius of upper cylinder
geo_yLower = 0.665; % height of lower conical part
geo_yPipeLower = 0.073; % lower height of oulet gas pipe
%% generate dInjection and traj_dInjection, generally should habe been done already
%% check input for calculating one given diameter
if isempty( find(dInjection == d) ) == 1
    error('Please enter a value belonging to injection diameter classes.');
end
%% Retrieve data for user-input diameter class
idx_d = find(dInjection == d);
data = traj_dInjection{idx_d,1};
%% Search for incomplete trajectories
dataSingleTraj = sortSingleTraj(data); %data for each trajectory
%% check if there's no incomplete trajectory in the given diameter class
checkIncomplete = zeros(1);
for i = 1:size(dataSingleTraj,1)
    checkIncomplete(i,1) = size(dataSingleTraj{i,1},1);
end
if isempty( find(checkIncomplete >= maxNumIter) ) == 1
    error('There are no incomplete trajectories in the given diameter class. Please try another diameter.');
end
%% calculate wall distance for incomplete trajectories
distanceWallArray = zeros(1); % initialize array for calulating distance to wall
for i = 1:size(dataSingleTraj,1)
    if size(dataSingleTraj{i,1},1) >= maxNumIter % find incomplete trajecotries
        distanceWallArray(i,1) = i - 1; % trajectory ID
        distanceWallArray(i,2:4) = dataSingleTraj{i,1}(end, 1:3); % endpoint x, y, z coordinates
        if distanceWallArray(i,3) > geo_yLower % check the height level of end point
            distanceWallArray(i,5) = geo_rUpper; % cylindrical part of chamber
        else
            distanceWallArray(i,5) = distanceWallArray(i,3) * geo_rUpper / geo_yLower; % conical part of chamber
        end    
        distanceWallArray(i,6) = distanceWallArray(i,5) - ... % end point distance to wall
            sqrt( distanceWallArray(i,2)^2 + distanceWallArray(i,4)^2 );
    end
end
distanceWall = distanceWallArray;
distanceWall(find(distanceWall(:,5)== 0),:) = []; % remove rows remained due to complete trajectories
%% print results
numInTraj = zeros(0);
numStream = 25;
fprintf('List of incomplete trajectories of injection diameter %g \x03bcm: \n', d*1e6)
for i = 1:size(distanceWall,1)
     if distanceWall(i,3) < geo_yPipeLower & distanceWall(:,6) < 0.05
         fprintf('Trajectory ID %g ends on height %.3g m. Probably escape as coarse product. \n', distanceWall(i,1), distanceWall(i,3) );
         numInTraj(i,1) = (idx_d - 1) * numStream + distanceWall(i,1) + 1;
     else
          fprintf('Trajectory ID %g ends on height %.3g m. Not sure. \n', distanceWall(i,1), distanceWall(i,3) );
     end
end
numInTraj = numInTraj(numInTraj ~= 0);
fprintf('Indices of "complete" trajectories in cell "traj": ['); 
fprintf(num2str(numInTraj(:).'));
fprintf('] \n \n');
end