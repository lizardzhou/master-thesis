function [dInjection, traj_dInjection] = sortDiamTraj(numStream, traj)
% Sort particle trajectories in cell traj according to diameters at
% injection. Sorted trajectories are stored in traj_dClass.
% Only applicable for same stream # for all injection diameters.
% Input:
%   numStream: # of streams per diameter
%   traj: cell containing info of all trajectories, generated by loadAllTrajectories.m. Call this function only after loadAllTrajectories.m is called. 
%   traj contains the following info:
%       (:,1) x
%       (:,2) y
%       (:,3) z
%       (:,4) residence time
%       (:,5) particle diameter calculated by Fluent (only used to check injection diameter)
%       ... (not relavant for calculating residence time and relative velocity)
%       (:,12) particle u_x
%       (:,13) particle u_y
%       (:,14) particle u_z
%       (:,15) gas u_x
%       (:,16) gas u_y
%       (:,17) gas u_z
%       (:,18) time step (cumulative sum of time step is residence time)
% Returns:
%   dInjection: lists all diameter values at injection
%   traj_dClass: cell of size numClass*1. Trajectories with same injection
%   diameter are stored in stacked form (trajectory ID (i+1) below
%   trajectory ID i) in each cell element.

% by Xiye Zhou, Oct. 2020
%% List the diameter values at injection
dInjection = zeros(size(traj,1)/numStream,1); % initialize array for diameters at injection
idxClass = 1;
for i = 1:numStream:size(traj,1)-1
    dInjection(idxClass,1) = traj{i,1}(1,5);
    idxClass = idxClass + 1;
end
%% Sort cell elements of same injection diameter into one cell element
traj_dInjection = cell(size(dInjection,1),1);
traj_idx = 1; % counts the position of first trajectory in an injection diameter
traj_idxLater = 2; % counts the index of traj which should be added to initialized traj_dClass{traj_idx,1}
for traj_diamIdx = 1:size(dInjection,1) % counts the index of injection diameters
    % initialize traj_dClass with first trajectory of every diameter
    traj_dInjection{traj_diamIdx,1} = traj{traj_idx,1};
    traj_idx = traj_idx + numStream;
    % add 2th to 25th trajectory below first trajectory for every injection diameter
    for i = traj_idxLater : (traj_idxLater+23)
         traj_dInjection{traj_diamIdx,1} = [traj_dInjection{traj_diamIdx,1};traj{i,1}]; 
     end
    traj_idxLater = traj_idxLater + numStream;
end
%% Retrieve injection diameter info
diam = num2str(dInjection' * 1e6);
fprintf('Injection diameter classes in \x03bcm: \n');
fprintf(diam);
fprintf('\nEach diameter class contains %g trajectories. \n \n', numStream);
end

