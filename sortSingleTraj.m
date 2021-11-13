function [traj] = sortSingleTraj(data)
% Extract single trajectories from array "data" containing all trajectories. Each single
% trajectory is stored in each element of cell "traj"
% Input:
%   data: cell containing all trajectories
% Return:
%   traj: stores single trajectories in each cell elements

% by Xiye Zhou, Oct. 2021
%%
% starting point of a trajectory is always 1.385 m high & in the center (y = 1.385, x = z = 0)
startPoint = find( (abs(data(:,2)-1.385)<1e-8) & (abs(data(:,1))<1e-5) & (abs(data(:,3))<1e-5) );
% every single trajectory is filled in the cell traj
traj = {data(startPoint(1):startPoint(2)-1,:)};
for i = 2:size(startPoint,1)-1
    traj(i,:) = {data(startPoint(i):startPoint(i+1)-1,:)};
end
traj(size(startPoint,1),:) = { data(startPoint(end):end,:) };
end

