% Lists indices of incomplete trajectories in cell "traj", which ends at
% bottom and can be treated as complete ones.
% by Xiye Zhou, Oct. 2021
% (:,1) injection diameter
% (:,2) # of incomplete trajectories that can be treated complete
% (:,3) index of incomplete trajectories in "traj"
%% Sort trajectories according to injection diameters
numStream = 25;
[dInjection, traj_dInjection] = sortDiamTraj(numStream, traj);
%% Listing
stat = cell(1,3);
for i = 8:16
    [~, ~, numInTraj] =calcPosIncompleteEnd(dInjection(i), 150e3, dInjection, traj_dInjection);
    stat{i-7,1} = dInjection(i);
    stat{i-7,2} = size(numInTraj,1);  
    stat{i-7,3} = numInTraj;
end