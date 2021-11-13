function [idxAcceptComplete, outDataPSDpAccept] = calcPosIncompleteEndAllDiam(maxNumIter, injectData, traj)
% Similar to function calcPosIncompleteEnd.m, calculate the distance to
% chamber wall of the end point of incomplete trajectories for all
% diameters in cell "traj" and get incomplete trajectories having end point
% at bottom, thus can be treated as coarse product in PSD calculation
% Input:
%   maxNumIter: max. # of data points per trajectory
%   injectData: diameter and mass flow of injections
%   traj: cell containing all trajectories, generated by
%   loadAllTrajectories.m
% Return:
%   outDataPSDpAccept: information of incomplete trajectories, which can be
%   treated as complete ones as coarse product, includes:
%       (:,1) y
%       (:,2) particle diameter at start point
%       (:,3) particle diameter at end point
%       (:,4) trajectory mass flow at start point
%       (:,5) single particle mass at start point
%       (:,6) single particle mass at end point
%       (:,7) trajectory mass flow at end point
%   idxAcceptComplete: # of such trajectories in cell "traj"

% by Xiye Zhou, Oct. 2021
%% Geometry of drying chamber
geo_rUpper = 0.397; % radius of upper cylinder
geo_yLower = 0.665; % height of lower conical part
geo_yPipeLower = 0.073; % lower height of oulet gas pipe
%% Sort out incomplete trajectories and calculate wall distance
distanceWallArrayAll = zeros(1); % initialize array for calulating distance to wall
for i = 1:size(traj,1)
    if size(traj{i,1},1) >= maxNumIter % find incomplete trajecotries
        distanceWallArrayAll(i,1) = i; % trajectory #
        distanceWallArrayAll(i,2:4) = traj{i,1}(end, 1:3); % endpoint x, y, z coordinates
        if distanceWallArrayAll(i,3) > geo_yLower  % check the height level of end point
            distanceWallArrayAll(i,5) = geo_rUpper; % cylindrical part of chamber
        else
            distanceWallArrayAll(i,5) = distanceWallArrayAll(i,3) * geo_rUpper / geo_yLower; % conical part of chamber
        end    
        distanceWallArrayAll(i,6) = distanceWallArrayAll(i,5) - ... % end point distance to wall
            sqrt( distanceWallArrayAll(i,2)^2 + distanceWallArrayAll(i,4)^2 );
    end
end
distanceWallArrayAll = distanceWallArrayAll( distanceWallArrayAll(:,1)~= 0,:); % remove rows remained due to complete trajectories
%% sort out trajectories ending at bottom and can be considered complete as coarse product
idxAcceptComplete = distanceWallArrayAll( find(distanceWallArrayAll(:,3) < geo_yPipeLower & ...
                                                                                    distanceWallArrayAll(:,6) < 0.05), : );
% if treat all incomplete trajectories as coarse product...
% idxAcceptComplete = distanceWallArrayAll(:,1);
%%
%%%%%%%%%%
%  Take a break   %
%%%%%%%%%%
%%  generate array for calculating PSD
outDataPSDpAccept = zeros(1);
for i = 1:size(idxAcceptComplete,1)
    outDataPSDpAccept(i,1) = traj{idxAcceptComplete(i),1}(end,2); % y
    outDataPSDpAccept(i,2) = traj{idxAcceptComplete(i),1}(1,6); % start diameter
    outDataPSDpAccept(i,3) = traj{idxAcceptComplete(i),1}(end,6); % end diameter
    for j = 1: size(injectData,1)
        if abs(outDataPSDpAccept(i,2) - injectData(j,1)) < 1e-8
            outDataPSDpAccept(i,4) = injectData(j,3); % start trajectory mass flow
        end
    end
    outDataPSDpAccept(i,5) = traj{idxAcceptComplete(i),1}(1,10); % start particle mass
    outDataPSDpAccept(i,6) = traj{idxAcceptComplete(i),1}(end,10); % end particle mass
    outDataPSDpAccept(i,7) = outDataPSDpAccept(i,6) * outDataPSDpAccept(i,4) / outDataPSDpAccept(i,5); % end trajectory mass flow
end
end