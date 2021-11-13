%% Script for calculating distribution of complete trajectories
% 1 - load and process data of injection
% 2 - calculate PSD of mass flow (q_3, Q_3), Sauter diameter, mode and median for injection
% 3 - calculate PSD of number (q_0, Q_0) for injection
% 4 - plot PSD for injection
% 5 - load and process data of all trajectories
% 6 - calculate mass flow-averaged residence time and relative velocity of
% all complete trajectories in whold spray dryer
% 7 - sort out first point "outStart" and last point "outEnd" of each complete trajectory
% 8 - calculate mass flow of each complete trajectory
% 9 - separate trajectories into fine and coarse product
% *calculate and plot residence time and relative velocity distribution (q_3, Q_3) for fine
% and coarse products
% 10 - calculate and plot PSD of mass flow (q_3, Q_3) for fine product
% 11 - calculate and plot PSD of mass flow (q_3, Q_3) for coarse product
% 12 - plot q and Q for injection and products

% by Xiye Zhou, Oct. 2021
%% Calculate and plot PSD of injection
%1 - load injection data
numStream = 25; % # of streams in each diameter class
rawDataInjection = readtable('0-injection.txt','PreserveVariableNames',true);
injectData = table2array(rawDataInjection);
% (:,1) particle diameter [micron]
% (:,2) mass flow rate of the whole diameter class [kg/s]
% (:,3) mass flow rate of each stream [kg/s]
injectData(:,1) = injectData(:,1) * 1e-6; % convert micron to meter
injectData(:,3) = injectData(:,2) / numStream; % mass flow of a single stream in each diameter class

% 2 - calculate PSD of mass flow for injection 
inCalcPSD(:,1) = [0;injectData(1:end-1,1)]; % lower limit of diameter
inCalcPSD(:,2) = injectData(:,1); % upper limit of diameter
inCalcPSD(:,3) = (inCalcPSD(:,1) + inCalcPSD(:,2)) / 2; % d_m,i
inCalcPSD(:,4) = inCalcPSD(:,2) - inCalcPSD(:,1); % Delta d_i
inFlowTotal = sum(injectData(:,2)); % total mass flow
inCalcPSD(:,5) = injectData(:,2); % massFlow_i
inCalcPSD(:,6) = inCalcPSD(:,5) ./ (inFlowTotal * inCalcPSD(:,4)); % q_3,i
inCalcPSD(:,7) = cumsum(inCalcPSD(:,5) / inFlowTotal); % Q_3,i
inCalcPSD(:,8) = inCalcPSD(:,3).^(-1) .* inCalcPSD(:,4) .* inCalcPSD(:,6); % summand of M_-1,3
in_d32 = 1 / sum(inCalcPSD(:,8)); % Sauter diameter = 1/(M_-1,3)
[~,in_idxMax] = max(inCalcPSD(:,6)); % index of the maximum q_3,i
in_dMode = inCalcPSD(in_idxMax,3); % mode diameter
in_dMedian = interp1(inCalcPSD(:,7),inCalcPSD(:,2),0.5); % median diameter
fprintf('Sauter diameter at injection is %.3g \x03bcm. \n',in_d32*1e6);
fprintf('Mode diameter at injection is %.3g \x03bcm. \n', in_dMode*1e6);
fprintf('Median diameter at injection is %.3g \x03bcm. \n \n', in_dMedian*1e6);

% 3 - calculate PSD of number for injection 
% Conversion from type 3 (mass flow) to type 0 (number): q_0,i = q_3,i * d_m,i^(-3) / M_-3,3
inCalcPSD(:,9) = inCalcPSD(:,3).^(-3) .* inCalcPSD(:,4) .* inCalcPSD(:,6); % summand of M_-3,3
inCalcPSD(:,10) = inCalcPSD(:,6) .* inCalcPSD(:,3).^(-3) / sum(inCalcPSD(:,9)); %q_0,i
inCalcPSD(:,11) = inCalcPSD(:,10) .* inCalcPSD(:,4); % summand of Q_0,i
inCalcPSD(:,12) = cumsum(inCalcPSD(:,11)); % Q_0,i
%% 4 - plot PSD of mass flow for injection
figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,3),inCalcPSD(:,6),'o-');
grid on
% title('Density distribution of mass flow at injection');
xlabel('Particle diameter $d_{m,i}$ [m]','Interpreter','latex');
ylabel('$q_3$ [m$^{-1}$]','Interpreter','latex');
figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,2),inCalcPSD(:,7),'o-');
grid on
% title('Cumulative distribution of mass flow at injection');
xlabel('Particle diameter $d_{max,i}$ [m]','Interpreter','latex');
ylabel('$Q_3$ [$-$]','Interpreter','latex');

% 4 - plot PSD of number for injection
figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,3),inCalcPSD(:,10),'o-');
grid on
% title('Density distribution of number at injection');
xlabel('Particle diameter $d_{m,i}$ [m]','Interpreter','latex');
ylabel('$q_0$ [m$^{-1}$]','Interpreter','latex');
% set(gca,'Xscale','log')
figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,2),inCalcPSD(:,12),'o-');
grid on
% title('Cumulative distribution of number at injection');
xlabel('Particle diameter $d_{max,i}$ [m]','Interpreter','latex');
ylabel('$Q_0$ [$-$]','Interpreter','latex');
% set(gca,'Xscale','log')
%%
%%%%%%%%%%
%  Take a break   %
%%%%%%%%%%
%% 5 - load data for all particle trajectories, generally should have been executed already
% loadAllTrajectories;
%% 6 - calculate mass flow-averaged values of all complete trajectories 
maxNumIter = 150e3;
[avgTauOut10k, avgRelVelOut10k, avgFracWater10k, outDataPSD, outStart, outEnd] = ...
    calcAvgOut(maxNumIter, traj, injectData); % max. # of data points from Fluent
% optional, for checking how to define "complete" to make residence time
% corresponds with experimental results
% [avgTauOut100k, avgRelVelOut100k, avgFracWater100k, outDataAvg100k,~,~] = calcAvgOut(100e3, traj, injectData);
% [avgTauOut50k, avgRelVelOut50k, avgFracWater50k, outDataAvg50k,~,~] = calcAvgOut(50e3, traj, injectData);
% [avgTauOut10k, avgRelVelOut10k, avgFracWater10k, outDataAvg10k,~,~] = calcAvgOut(10e3, traj, injectData);
%% 7 - sort out first point "outStart" and last point "outEnd" of each complete trajectory
% outStart and outEnd are obtained from step 6
%% 8 - generate array "outDataPSD" for calculating PSD of mass flow for each complete trajectory
% outDataPSD is obtained from step 6
%%
%%%%%%%%%%
%  Take a break   %
%%%%%%%%%%
%% 9 - sort out gas and particle outlets from array outDataPSD into outDataPSDg and outDataPSDp
outCountP = 1;
outCountG = 1;
outDataPSDp = zeros(size(outEnd,1),size(outDataPSD,2)); % initialize array for escaped coarse particles
outDataPSDg = zeros(size(outEnd,1),size(outDataPSD,2)); % initialize array for escaped fine particles
for i = 1:size(outEnd,1)
    if outDataPSD(i,1) < 1e-5 % coarse product always escapes from bottom of drying chamber (y = 0)
        outDataPSDp(outCountP,:) = outDataPSD(i,:);
        outCountP = outCountP + 1;
    else
        outDataPSDg(outCountG,:) = outDataPSD(i,:);
        outCountG = outCountG + 1;
    end
    outDataPSDp = outDataPSDp(1:outCountP-1,:);
    outDataPSDg = outDataPSDg(1:outCountG-1,:);
end
% find incomplete trajectories ending at bottom part of drying chamber,
% they can be considered complete as coarse product
[idxAcceptComplete, outDataPSDpAccept] = calcPosIncompleteEndAllDiam(maxNumIter, injectData, traj);
% print out results
fprintf('%g trajectories escaped as fine product. \n',size(outDataPSDg,1));
fprintf('%g trajectories escaped as coarse product (including %g incomplete trajectories ending near particle outlet). \n \n', ...
    size(outDataPSDp,1)+size(idxAcceptComplete,1), size(idxAcceptComplete,1));
%% *calculate Q_3 according to residence time and relative velocity
% for fine product, tau50
outDataPSDgTau(:,1) = outDataPSDg(:,8); % residence time at end point
outDataPSDgTau(:,2) = outDataPSDg(:,7); % mass flow at end point
outDataPSDgTau = sortrows(outDataPSDgTau,1);
mFlowG = sum(outDataPSDgTau(:,2));
outDataPSDgTau(:,3) =  cumsum(outDataPSDgTau(:,2) / mFlowG); % Q_3
tau50g = interp1(outDataPSDgTau(:,3),outDataPSDgTau(:,1),0.5);
fprintf('It takes %.3g s to get half of mass flow of fine product out from gas outlet. \n', tau50g);
% for fine product, relVel50
outDataPSDgRelVel(:,1) = outDataPSDg(:,10); % relative velocity at end point
outDataPSDgRelVel(:,2) = outDataPSDg(:,7); % mass flow at end point
outDataPSDgRelVel = sortrows(outDataPSDgRelVel,1);
outDataPSDgRelVel(:,3) =  cumsum(outDataPSDgRelVel(:,2) / mFlowG); % Q_3
relVel50g = interp1(outDataPSDgRelVel(:,3),outDataPSDgRelVel(:,1),0.5);
fprintf('Half of mass flow out from gas outlet has relative velocity lower than %.3g m/s. \n \n', relVel50g);

figure % plot tau
plot(outDataPSDgTau(:,1), outDataPSDgTau(:,3), 'o-');
set(gca,'Xscale','log');
% title('Distribution of mass flow according to residence time at gas outlet');
xlabel('Residence time [s]','Interpreter','latex');
ylabel('$Q_3$ [$-$]','Interpreter','latex');
ylim([0 1]);
grid on
figure % plot relVel
plot(outDataPSDgRelVel(:,1), outDataPSDgRelVel(:,3), 'o-');
% set(gca,'Xscale','log');
% title('Distribution of mass flow according to relative velocity at gas outlet');
xlabel('Relative velocity [m/s]','Interpreter','latex');
ylabel('$Q_3$ [$-$]','Interpreter','latex');
ylim([0 1]);
grid on

% for coarse product, tau50
outDataPSDpTau(:,1) = outDataPSDp(:,8); % residence time
outDataPSDpTau(:,2) = outDataPSDp(:,7); % mass flow at end point
outDataPSDpTau = sortrows(outDataPSDpTau,1);
mFlowP = sum(outDataPSDpTau(:,2));
outDataPSDpTau(:,3) =  cumsum(outDataPSDpTau(:,2) / mFlowP); % Q_3
tau50p = interp1(outDataPSDpTau(:,3),outDataPSDpTau(:,1),0.5);
fprintf('It takes %.3g s to get half of mass flow of coarse product out from particle outlet. \n', tau50p);
% for coarse product, relVel50
outDataPSDpRelVel(:,1) = outDataPSDp(:,10); % relative velocity at end point
outDataPSDpRelVel(:,2) = outDataPSDp(:,7); % mass flow at end point
outDataPSDpRelVel = sortrows(outDataPSDpRelVel,1);
outDataPSDpRelVel(:,3) =  cumsum(outDataPSDpRelVel(:,2) / mFlowP); % Q_3
relVel50p = interp1(outDataPSDpRelVel(:,3),outDataPSDpRelVel(:,1),0.5);
fprintf('Half of mass flow out from particle outlet has relative velocity lower than %.3g m/s. \n \n', relVel50p);

figure % plot tau
plot(outDataPSDpTau(:,1), outDataPSDpTau(:,3),'o-');
% title('Distribution of mass flow according to residence time at particle outlet');
xlabel('Residence time [s]','Interpreter','latex');
ylabel('$Q_3$ [$-$]','Interpreter','latex');
ylim([0 1]);
grid on
figure % plot relVel
plot(outDataPSDpRelVel(:,1), outDataPSDpRelVel(:,3), 'o-');
% title('Distribution of mass flow according to relative velocity at particle outlet');
xlabel('Relative velocity [m/s]','Interpreter','latex');
ylabel('$Q_3$ [$-$]','Interpreter','latex');
ylim([0 1]);
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10 - calculate PSD of fine product using array outCalcPSDg
diamSortedG = sortrows(outDataPSDg(:,3)); % sort diameter ascendingly in order to see the diameter class
listClassG = diamSortedG * 1e6; % list overview of outlet diameters
% take a look at listClassG to adjust diameter classes for calcPSDout!!!
[outG_d32, outG_dMode, outG_dMedian, outCalcPSDg] = ...
    calcPSDout( [0,1:2:15,19:4:35]*1e-6, 41e-6, outDataPSDg, 'fine');
%% 11 - calculate PSD of coarse product using array outCalcPSDpAll
% data for complete trajectories
outDataPSDpAll = outDataPSDp(:,1:7);
if isempty(idxAcceptComplete) == 1 % check if there are acceptable incomplete trajectories
    outDataPSDpAll = outDataPSDp; 
else
    outDataPSDpAll = [outDataPSDpAll; outDataPSDpAccept]; % combine data for incomplete but acceptable trajectories
end

% calculate PSD 
diamSortedP = sortrows(outDataPSDpAll(:,3)); % sort diameter ascendingly in order to see the diameter class
listClassP = diamSortedP * 1e6; % list overview of outlet diameters
% take a look at listClassG to adjust diameter classes for calcPSDout!!!
[outP_d32, outP_dMode, outP_dMedian, outCalcPSDpAll] = ...
    calcPSDout( [0,20:3:38]*1e-6, 41e-6, outDataPSDpAll, 'coarse');
%% 12 - plot PSD for injection, fine and coarse products
figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,3),inCalcPSD(:,6),'o-', ...
    outCalcPSDg(:,3),outCalcPSDg(:,6),'o-', ...
    outCalcPSDpAll(:,3),outCalcPSDpAll(:,6),'o-');
grid on;
% title('Density distribution of mass flow');
xlabel('$d_{m,i}$ [m]','Interpreter','latex');
ylabel('$q_3$ [m$^{-1}$]','Interpreter','latex');
legend('Injection','Fine product','Coarse product','Interpreter','latex','Location','best');
figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,2),inCalcPSD(:,7),'o-', ...
    outCalcPSDg(:,2),outCalcPSDg(:,8),'o-', ...
    outCalcPSDpAll(:,2),outCalcPSDpAll(:,8),'o-');
grid on;
% title('Cumulative distribution of mass flow');
xlabel('$d_{max,i}$ [m]','Interpreter','latex');
ylabel('$Q_3$ [$-$]','Interpreter','latex');
legend('Injection','Fine product','Coarse product','Interpreter','latex','Location','best');

figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,3),inCalcPSD(:,10),'o-', ...
    outCalcPSDg(:,3),outCalcPSDg(:,11),'o-', ...
    outCalcPSDpAll(:,3),outCalcPSDpAll(:,11),'o-');
set(gca,'Xscale','log');
grid on;
% title('Density distribution of number');
xlabel('$d_{m,i}$ [m]','Interpreter','latex');
ylabel('$q_0$ [m$^{-1}$]','Interpreter','latex');
legend('Injection','Fine product','Coarse product','Interpreter','latex','Location','best');
figure
set(gcf,'renderer','Painters')
plot(inCalcPSD(:,2),inCalcPSD(:,12),'o-', ...
    outCalcPSDg(:,2),outCalcPSDg(:,13),'o-', ...
    outCalcPSDpAll(:,2),outCalcPSDpAll(:,13),'o-');
set(gca,'Xscale','log');
grid on;
% title('Cumulative distribution of number');
xlabel('$d_{max,i}$ [m]','Interpreter','latex');
ylabel('$Q_0$ [$-$]','Interpreter','latex');
legend('Injection','Fine product','Coarse product','Interpreter','latex','Location','best');