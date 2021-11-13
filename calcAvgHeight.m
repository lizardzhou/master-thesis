function [tauAvg, tauDev, tauArray, relVelAvg, relVelArray] = calcAvgHeight(d, maxNumIter, dInjection, traj_dInjection, hLow, hHigh, plotAll)
% Calculate average residence time and relative velocity between two height levels of a given
% diameter class in the bulk space of chamber
% Bulk space: space above gas outlet pipe, 
% Input:
%   d: user-input injection diameter in meter
%   maxNumIter: maximal number of data points per trajectory
%   dInjection: lists all injection diameters ascendingly. Generated in
%   dataAvg.m by sortDiamTraj.m
%   traj_dInjection: cell of size numClass*1. Each cell element contains
%   particle of same injection diameter. Generated in dataAvg.m by sortDiamTraj.m
%   hLow: lower height in meter
%   hHigh: upper height in meter
%   plotAll: if user sets the value as 1, parameters will be plotted; in case of other
%   values, no plot will be made
% Return:
%   tauAvg: arithmetic average of residence time difference between hLow and hHigh
%   tauDev: standard deviation of residence time difference between hLow and hHigh
%   tauArray: calculation info of residence time difference of each trajectory between hLow
%   and hHigh
%   relVelAvg: average particle relative velocity to gas calculated from
%   time integral
%   relVelArray: calculation info of relative velocity 

% by Xiye Zhou, Oct. 2021
%% Geometry of drying chamber
hGasOutHigh = 0.125; % highest height of gas outlet pipe
%% Check input
% Check height
if hLow < 0 || hHigh < 0
    error('Please enter non-nagative numbers!');
else
    if hLow > 1.52 || hHigh > 1.52
        error('Height value should not be higher than the tower height 1.52 m!');
    end
end
if hHigh == hLow
    error('Please enter two different height values.');
else
    if hHigh < hLow % swap hHigh and hLow values if user enters them inversely
        exchg = hHigh;
        hHigh = hLow;
        hLow = exchg;
    end
end
% check d
if isempty( find(dInjection == d) ) == 1
    error('Please enter a value belonging to injection diameter classes.');
end
% in case of outlet regions
if hLow < hGasOutHigh && hHigh > hGasOutHigh
    warning('Height values less than 0.126 m belong to outlet regions and is calculated in dataAvg.m. Self-defined interval should be above 0.126 m. Otherwise the result may become inaccurate due to outflows.');
end
%% Retrieve data for user-input diameter class
idx_d = find(dInjection == d);
data = traj_dInjection{idx_d,1};
%% Calculate residence time for complete trajectories
dataSingleTraj = sortSingleTraj(data); %data for each trajectory
tauArray = zeros(size(dataSingleTraj,1),1); % initialize tauArray for calculating residence time for each trajectory
% calculate residence time for each trajectory
for i = 1:size(dataSingleTraj,1)
    tauArray(i,1) = i - 1; % trajectory ID in Fluent
    idxHeightTau = find( (dataSingleTraj{i,1}(:,2)>=hLow) & (dataSingleTraj{i,1}(:,2)<=hHigh) ); % indices of height between hLow and hHigh
    tauArray(i,2) = sum(dataSingleTraj{i,1}(idxHeightTau,18)); % sum of time steps between hLow and hHigh
end
% calculate average residence time between hLow and hHigh for complete
% trajectories
tauArray(:,3) = tauArray(:,2);
idxCompleteSpace = zeros(1); % count index for complete trajectories
idxIncompleteSpace = zeros(1); % count index for incomplete trajectories
for i = 1:size(tauArray,1)
    if size(dataSingleTraj{i,1},1) < maxNumIter
        idxCompleteSpace(i,1) = i;
    else
        idxIncompleteSpace(i,1) = i;
        tauArray(i,3) = -1;
    end
end
idxIncomplete = find(tauArray(:,3) == -1); % indices for incomplete trajectory without 0 spacing in idxIncomplete
idxComplete = find(tauArray(:,3) ~= -1); % indices for complete trajectory without 0 spacing in idxComplete
idxCompleteNonZero = find( (tauArray(:,3) ~= -1) & (tauArray(:,3) ~= 0) );
tauAvg = mean( tauArray(idxCompleteNonZero,3) );
tauDev = std( tauArray(idxCompleteNonZero,3) );

% plot and display residence time of each trajectory if the user wants to do so
if plotAll == 1
    % plot
    figure
    x0=500;
    y0=500;
    width=900;
    height=400;
    set(gcf,'position', [x0,y0,width,height])
    bCompleteTau = bar( tauArray(idxComplete,1), tauArray(idxComplete,2), 0.5  ); % comlete trajectories are marked blue
    bCompleteTau.FaceColor = '#0072BD'; % Hexadecimal Color Code for blue
    xlim([-1 25]);
    xticks(0:24);
    hold on
    bInompleteTau = bar( tauArray(idxIncomplete,1), tauArray(idxIncomplete,2), 0.5 ); % incomplete trajectories are marked red
    bInompleteTau.FaceColor = '#A2142F';  % Hexadecimal Color Code for red 
    title("Residence time for all trajectories of " + dInjection(idx_d)*1e6 + " $\mu$m between heights " +  hLow + " and " + hHigh + " m",'Interpreter','latex');
    xlabel('Trajectory ID in Fluent','Interpreter','latex');
    ylabel('Residence time [s]','Interpreter','latex');
    set(gca,'Yscale','log');
    xlim([-1 25]);
    xticks(0:24);
    ylim([1e-3 1e2]);
    hold off
    % display residence time result
    fprintf('%g of %g trajectories are complete and applied to calculate average residence time and relative velocity. \n', ...
    size(idxComplete,1), size(tauArray,1));
    fprintf('residence time between %g m and %g m of injection diameter %g \x03bcm: %.3g +- %.3g s \n', ...
    hLow, hHigh, dInjection(idx_d)*1e6, tauAvg, tauDev);
end
%% Calculate velocity relative of particle to gas for all trajectories
relVelArray = zeros(size(dataSingleTraj,1),1); % initialize relVelArray for calculating relative velocity for each trajectory
% Calculate relative velocity for each trajectory
for i = 1:size(dataSingleTraj,1)
    relVelArray(i,1) = i - 1; % trajectory ID in Fluent
    if tauArray(i,2) == 0 % if residence time is zero (outside the region), the velocity should also be zero
        relVelArray(i,2) = 0;
    else
        dataSingleTraj{i,1}(:,19) = sqrt( (dataSingleTraj{i,1}(:,12) - dataSingleTraj{i,1}(:,15)).^2 + ... % relative velocity of particle to gas
                                                              (dataSingleTraj{i,1}(:,13) - dataSingleTraj{i,1}(:,16)).^2 + ...
                                                              (dataSingleTraj{i,1}(:,14) - dataSingleTraj{i,1}(:,17)).^2 ); 
        idxHeightVel = find( (dataSingleTraj{i,1}(:,2)>=hLow) & (dataSingleTraj{i,1}(:,2)<=hHigh) ); % indices of height between hLow and hHigh
        relVelTimeStep = dataSingleTraj{i,1}(idxHeightVel,18); % time step values between hLow and hHigh
        relVelStep = dataSingleTraj{i,1}(idxHeightVel,19); % relative velocity values between hLow and hHigh
        sumRelVel = relVelStep(1:end-1) + relVelStep(2:end); % summand of relVel(n) + relVel(n+1)
        % integral of relative velocity over time step using trapezoidal rule
        sumIntegral = sum(0.5 * relVelTimeStep(2:end) .* sumRelVel); % integral: sum of trapezoid area of each time step
        sumTime = sum(relVelTimeStep); % sum of time steps
        relVelSingle = sumIntegral / sumTime; % average of relative velocity
        relVelArray(i,2) = relVelSingle;
    end
end
% calculate average residence time between hLow and hHigh for complete
% trajectories
relVelArray(:,3) = relVelArray(:,2);
for i = 1:size(relVelArray,1)
    if size(dataSingleTraj{i,1},1) > maxNumIter
        relVelArray(i,3) = -1;
     end
end
relVelAvg = mean( relVelArray(idxCompleteNonZero,3) );
relVelDev = std( relVelArray(idxCompleteNonZero,3) );

% plot and display relative velocity of each trajectory if the user wants to do so
if plotAll == 1
    % plot
    figure
    x0=500;
    y0=500;
    width=900;
    height=400;
    set(gcf,'position', [x0,y0,width,height])
    bCompleteRelVel = bar( relVelArray(idxComplete,1), relVelArray(idxComplete,2), 0.5 ); % complete trajectories are marked blue
    bCompleteRelVel.FaceColor = '#0072BD'; % Hexadecimal Color Code for blue
    xlim([-1 25])
    xticks(0:24);
    hold on
    bInompleteRelVel = bar( relVelArray(idxIncomplete,1), relVelArray(idxIncomplete,2), 0.5 ); % incomplete trajectories are marked red
    bInompleteRelVel.FaceColor = '#A2142F';  % Hexadecimal Color Code for red 
    title("Relative velocity for all trajectories of " + dInjection(idx_d)*1e6 + " $\mu$m between heights " +  hLow + " and " + hHigh + " m",'Interpreter','latex');
    xlabel('Trajectory ID in Fluent','Interpreter','latex');
    ylabel('Particle relative velocity to gas [m/s]','Interpreter','latex');
    set(gca,'Yscale','log');
    xlim([-1 25])
    xticks(0:24);
    ylim([1e-3 1e2])
    hold off   
    % display relative velocity result
    fprintf('relative velocity between %g m and %g m of injection diameter %g \x03bcm: %.3g +- %.3g m/s \n \n', ...
    hLow, hHigh, dInjection(idx_d)*1e6, relVelAvg, relVelDev);
end
end