%% Calculates the drying process of single trajectories
% by Xiye Zhou, Oct. 2021
%% Comparison of diameter of a single particle
[~, stau, sdFluent, sdUDF, smoi, sfracWFluent, sfracWUDF, sf, stemp, sdryingTime1, sdryingTime2, sdryingTime3] = ...
 calcDryingSingle(traj, 1, 25, 1);
%% temperature, moisture content, relative drying rate
% temperature
figure
yyaxis left
sptemp = plot(stau, stemp-273,'LineWidth', 1);
ylabel('Particle temperature [C$^\circ$]','Interpreter','latex');
hold on
yyaxis right
spmoi = plot(stau, smoi,'LineWidth',1);
ylabel('Moisture content [$-$]','Interpreter','latex');
xlabel('Residence time [s]','Interpreter','latex');
yline(0.7,':','$x_{cr}$','LineWidth',1.5,'Interpreter','latex');
yline(0.05,':','$x_{eq}$','LineWidth',1.5,'Interpreter','latex');
xline(0.08,':', 'LineWidth',1.5,'Interpreter','latex');
xline(0.12,':', 'LineWidth',1.5,'Interpreter','latex');
grid on
% Relative drying rate
figure
plot(stau, sf, 'LineWidth',1, 'Color', '#7E2F8E');
ylabel('Relative drying rate [$-$]','Interpreter','latex');
xlabel('Residence time [s]','Interpreter','latex');
xline(0.08,':', 'LineWidth',1.5,'Interpreter','latex');
xline(0.12,':', 'LineWidth',1.5,'Interpreter','latex');
grid on
%% Particle mass
figure
yyaxis left
plot(stau, sdFluent, stau, sdUDF,'LineWidth',1);
ylabel('Particle diameter [kg]','Interpreter','latex');
ylim([0 28e-6]);
hold on
yyaxis right
plot(stau, traj{381,1}(:,10),'LineWidth',1)
ylabel('Single particle mass [kg]','Interpreter','latex');
grid on
xlabel('Residence time [s]','Interpreter','latex');
legend('Diameter from Fluent','Diameter from UDF','Single particle mass','Interpreter','latex');
%% ID of trajectories to be studied
plotTrajID = [1, ... % d = 1
                      26, ... % d = 4
                      60, ... % d = 7
                      80, ... % d = 10
                      105, ... % d = 13
                      130, ... % d = 16
                      160, ... % d = 19
                      181, ... % d = 22
                      202, ... % d = 27
                      234, ... % d = 32
                      255, ... % d = 37
                      284, ... % d = 42
                      318, ... % d = 47
                      334, ... % d = 52
                      373, ... % d = 57
                      381]'; % d = 60
%% Calculate injection diameter, tau, temperature and f for each trajectory
dStream = zeros(size(plotTrajID)); %array for injection diameter
dTau = cell(size(plotTrajID)); % array for residence time
dTemp = cell(size(plotTrajID)); %array for particle temperature
df = cell(size(plotTrajID)); %array for relative drying rate
dDryingTime1 = zeros(size(plotTrajID));%duration of 1st drying stage
dDryingTime2 = zeros(size(plotTrajID)); %duration of 2nd drying stage
dDryingTime3 = zeros(size(plotTrajID)); %duration of 3rd drying stage
for i = 1:size(plotTrajID,1)
    [dStream(i,1),~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = calcDryingSingle(traj, plotTrajID(i), 25,0);
    [~,temporalTau, ~, ~, ~, ~, ~, ~, ~, ~, ~] = calcDryingSingle(traj, plotTrajID(i), 25,0);
    dTau{i,1} = temporalTau;
    [~, ~, ~, ~, ~, ~, ~, ~, temporalTemp,~,~] = calcDryingSingle(traj, plotTrajID(i), 25,0);
    dTemp{i,1} = temporalTemp;
    [~,~, ~, ~, ~, ~, ~, temporalF,~, ~,~] = calcDryingSingle(traj, plotTrajID(i), 25,0);
    df{i,1} = temporalF;
    [~,~, ~, ~, ~, ~, ~, ~, ~, dDryingTime1(i,1), ~,~] = calcDryingSingle(traj, plotTrajID(i), 25,0);
    [~,~, ~, ~, ~, ~, ~, ~, ~, ~, dDryingTime2(i,1),~] = calcDryingSingle(traj, plotTrajID(i), 25,0);
    [~,~, ~, ~, ~, ~, ~, ~, ~, ~,~, dDryingTime3(i,1)] = calcDryingSingle(traj, plotTrajID(i), 25,0);
end
%% Plot trajectories of given IDs
plotTraj(dStream, traj, plotTrajID);
lgdStream = dStream * 1e6; % diameter in microns for plot legend
%% Plot drying times of given IDs
figure
set(gcf,'renderer','Painters')
dDryingTime = [dDryingTime1, dDryingTime2, dDryingTime3];
aDrying = area(lgdStream, dDryingTime);
aDrying(1).FaceColor = '#EDB120';
aDrying(2).FaceColor = '#D95319';
aDrying(3).FaceColor = '#A2142F';
xlabel('Injection diameter [$\mu$m]','Interpreter','latex');
ylabel('Drying time [s]','Interpreter','latex');
legend('1st drying stage','2nd drying stage','3rd drying stage', 'Interpreter','latex','Location','best')
grid on
%% Comparison of all trajectories
figure
set(gcf,'renderer','Painters')
% set plot size
plt_x0=500;
plt_y0=500;
plt_width=900;
plt_height=400;
% plot temperature change during drying
% set plot colormap
colormap parula;
cmap = colormap;
colorCount = 1;
set(gcf,'position', [plt_x0,plt_y0,plt_width,plt_height])
for i = 1:size(plotTrajID,1)
    plotColor = cmap(colorCount, :);
    plot(dTau{i,1}, dTemp{i,1}-273, 'Color', plotColor,'LineWidth',1);
    hold on;
    colorCount = colorCount + size(cmap,1) / size(plotTrajID,1);
end
xlabel('Residence time [s]','Interpreter','latex');
xlim([0 1.05])
ylabel('Temperature [$^\circ$C]','Interpreter','latex');
lgdDiam = char(ones(size(lgdStream)) * '$$d_{inj}$$ = ');
lgdNum = num2str(lgdStream, '%d');
lgdUnit = char(ones(size(lgdStream)) * ' $$\mu$$m');
legend([lgdDiam lgdNum lgdUnit],'Interpreter','latex','Location','best');
grid on
% plot relative drying rate change during drying
figure
set(gcf,'renderer','Painters')
colorCount = 1;
set(gcf,'position', [plt_x0,plt_y0,plt_width,plt_height])
for i = 1:size(plotTrajID,1)
    plotColor = cmap(colorCount, :);
    plot(dTau{i,1}, df{i,1}, 'Color', plotColor,'LineWidth',1);
    hold on;
    colorCount = colorCount + size(cmap,1) / size(plotTrajID,1);
end
xlabel('Residence time [s]','Interpreter','latex');
ylabel('Relative drying rate [$-$]','Interpreter','latex');
lgdDiam = char(ones(size(lgdStream)) * '$$d_{inj}$$ = ');
lgdNum = num2str(lgdStream, '%d');
lgdUnit = char(ones(size(lgdStream)) * ' $$\mu$$m');
legend([lgdDiam lgdNum lgdUnit],'Interpreter','latex','Location','best');
grid on
%% Sort trajectories according to injection diameters
numStream = 25;
[dInjection, traj_dInjection] = sortDiamTraj(numStream, traj);
%% Calculate the mass-flow average drying times of all injection diameters
% Calculate average drying times for each injection diameter
allDryingTime = zeros(size(dInjection,1),4);
for i = 1: size(dInjection,1)
    [allDryingTime(i,1), ~, ~] = calcDryingTimeDiam(dInjection(i), 25, traj, dInjection);
    [~, allDryingTime(i,2), ~, ~, ~] = calcDryingTimeDiam(dInjection(i), 25, traj, dInjection);
    [~, ~, allDryingTime(i,3), ~,  ~] = calcDryingTimeDiam(dInjection(i), 25, traj, dInjection);
    [~, ~, ~, allDryingTime(i,4), ~] = calcDryingTimeDiam(dInjection(i), 25, traj, dInjection);
end
% Load injection data to average values by mass flow in array injectData
%   (:,1) particle diameter [micron]
%   (:,2) mass flow rate of the whole diameter class [kg/s]
rawDataInjection = readtable('0-injection.txt','PreserveVariableNames',true);
injectData = table2array(rawDataInjection);
flowInject = sum(injectData(:,2));
% Calculate mass-flow average of the drying times
avgDryingTime1 = sum( allDryingTime(:,1) .* injectData(:,2) ./ flowInject);
avgDryingTime2 = sum( allDryingTime(:,2) .* injectData(:,2) ./ flowInject);
%% Plot
figure
set(gcf,'renderer','Painters')
avgDryingTime = [allDryingTime(:,1), allDryingTime(:,3) ];
avgDrying = area(dInjection, avgDryingTime);
avgDrying(1).FaceColor = '#EDB120';
avgDrying(2).FaceColor = '#D95319';
xlabel('Injection diameter [m]','Interpreter','latex');
ylabel('Average drying time [s]','Interpreter','latex');
legend('1st drying stage','2nd drying stage', 'Interpreter','latex','Location','best')
grid on
% Plot for fitting
figure
fit1 = fit(injectData(:,1), allDryingTime(:,1), 'poly4','normalize','on');
plot(fit1, injectData(:,1), allDryingTime(:,1))
hold on 
fit2 = fit(injectData(:,1), allDryingTime(:,3), 'poly4','normalize','on');
plot(fit2, injectData(:,1), allDryingTime(:,3))
hold on
errorbar(injectData(:,1), allDryingTime(:,1), allDryingTime(:,2),'o', 'Color', '#0072BD');
hold on
errorbar(injectData(:,1), allDryingTime(:,3), allDryingTime(:,4),'*', 'Color', '#0072BD');
xlabel('Injection diameter [$\mu$m]','Interpreter','latex');
ylabel('Average drying time [s]','Interpreter','latex');
ylim([0 0.6]);
legend('','Fitting','','','1st drying stage','2nd drying stage', 'Interpreter','latex','Location','best')
grid on