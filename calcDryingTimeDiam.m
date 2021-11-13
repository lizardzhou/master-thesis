function [allTimeAvg1, allTimeDev1, allTimeAvg2, allTimeDev2, allDryingTime] = calcDryingTimeDiam(d, numStream, traj, dInjection)
% Calculate the 1st and 2nd drying times for all trajectories with injection diameter d
% Input:
%   d: injection diameter given in meter
%   numStream: # of streams in each injection diameter
%   traj: cell "traj" genenrated by loadAllTrajctories.m
%   dInjection: list of injection diameters
% Return:
%   allTimeAvg1: average drying time of 1st drying stage
%   allTimeDev1: standard deviation of allAvgTime1
%   allTimeAvg2: average drying time of 2nd drying stage
%   allTimeDev2: standard deviation of allAvgTime2
%   allDryingTime: drying times of all trajectories
%       (:,1) time of 1st drying stage
%       (:,2) time of 2nd drying stage

% by Xiye Zhou, Nov. 2021
%% Calculate the average drying time for 1st and 2nd stage
allDryingTime = zeros(numStream, 2);
idx_d = find(dInjection == d);
idx_dTraj = (idx_d - 1) * 25 + 1;
for i = idx_dTraj:idx_dTraj+24
    singleData = traj{i,1};
    f = singleData(:,9); % relative drying rate
    idx_f1 = find(f<1,1); % begin of 2nd drying stage: first f-value smaller than 1
    idx_f0 = find(f<5e-5,1); % end of 2nd drying stage: first f-value smaller than 5e-5 (approach 0)
    allDryingTime(i-(idx_d-1)*25,1) = sum(singleData(1:idx_f1,18)); % duration of 1st drying stage
    allDryingTime(i-(idx_d-1)*25,2) = sum(singleData(idx_f1:idx_f0,18)); % duration of 2nd drying stage
end
allTimeAvg1 = mean( allDryingTime(:,1) );
allTimeDev1 = std(allDryingTime(:,1));
allTimeAvg2 = mean( allDryingTime(:,2) );
allTimeDev2 = std(allDryingTime(:,2));
end