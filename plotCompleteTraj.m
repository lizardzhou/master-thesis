maxNumIter = 150e3;
%% Sort trajectories according to injection diameters
numStream = 25;
[dInjection, traj_dInjection] = sortDiamTraj(numStream, traj);
%% Plot all complete trajectories
figure
% set plot size
plt_x0 = 300;
plt_y0 = 300;
plt_width = 350;
plt_height = 600;
set(gcf,'position', [plt_x0,plt_y0,plt_width,plt_height])
% plot all complete trajectories
for i = 1:size(traj,1)
    if size(traj{i,1}) < maxNumIter
        plot3(traj{i,1}(:,1),traj{i,1}(:,3),traj{i,1}(:,2) ,'Color','#0072BD');    % color option: ,'Color','#EDB120'
    end
    hold on;
end
grid on;
xlim([-0.397 0.397])
ylim([-0.397 0.397])
zlim([0 1.52])
xlabel('$x$ [m]', 'Interpreter','latex');
ylabel('$z$ [m]', 'Interpreter','latex');
zlabel('$y$ [m]', 'Interpreter','latex');
%% Plot trajectories colored with relative velocity / Stokes number for a single trajectory
viscosityG = 1.72e-5;
densityG = 1.039;
idxTraj = 299; % 279-fine, 30s  % 299-coarse, 69s
z_d = 42e-6;
zVelG = sqrt(  traj{idxTraj,1}(:,15).^2 + traj{idxTraj,1}(:,16).^2 + traj{idxTraj,1}(:,17).^2 );
zRelVel = sqrt( ( traj{idxTraj,1}(:,12) -  traj{idxTraj,1}(:,15) ).^2 + ...
                          ( traj{idxTraj,1}(:,13) -  traj{idxTraj,1}(:,16) ).^2 + ...
                          ( traj{idxTraj,1}(:,14) -  traj{idxTraj,1}(:,17) ).^2  );
zReRel = z_d .* zRelVel * densityG / viscosityG; 
zDragCoeff = zeros(size(traj{idxTraj,1},1),1);
for zi = 1:size(traj{idxTraj,1},1)
    if zReRel(zi) <0.1
        zDragCoeff(zi) = 24 / zReRel(zi);
    else
        if zReRel(zi) < 1
            zDragCoeff(zi) = 22.73 ./ zReRel(zi) + 0.0903 ./ zReRel(zi).^2 + 3.69;
        else
            if zReRel(zi) < 10
                zDragCoeff(zi) = 29.1667 ./ zReRel(zi) - 3.8889 ./ zReRel(zi).^2 + 1.222;
            else
                if zReRel(zi) < 100
                    zDragCoeff(zi) = 46.5 ./ zReRel(zi) - 116.67 ./ zReRel(zi).^2 + 0.6167;
                else
                    zDragCoeff(zi) = 98.33 ./ zReRel(zi) - 2778 ./ zReRel(zi).^2 + 0.3644;
                end
            end
        end
    end
end
zRelaxTau = densityG .* z_d.^2 * 24 ./ (18 * viscosityG .* zDragCoeff .* zReRel);
zStk = zRelaxTau .* zVelG ./ z_d;
figure
plt_widthRelVel = 500;
set(gcf,'position', [plt_x0,plt_y0,plt_widthRelVel,plt_height])
set(gcf,'renderer','Painters')
scatter3(traj{idxTraj,1}(:,1),traj{idxTraj,1}(:,3),traj{idxTraj,1}(:,2), 1, zStk)
grid on
xlim([-0.397 0.397])
ylim([-0.397 0.397])
zlim([0 1.52])
xlabel('$x$ [m]', 'Interpreter','latex');
ylabel('$z$ [m]', 'Interpreter','latex');
zlabel('$y$ [m]', 'Interpreter','latex');
zCBar = colorbar;
colormap turbo
zCBar.Label.String = 'Stokes number [$-$]';
zCBar.Label.Interpreter = 'latex';
caxis([0 1.8]);
%% plot complete trajectories of a given injection diameter
z_d = 19e-6;
idx_d = find(dInjection == z_d);
idx_dTraj = (idx_d - 1) * 25 + 1;
figure
% set plot size
set(gcf,'position', [plt_x0,plt_y0,plt_width,plt_height]);
for i = idx_dTraj:idx_dTraj+24
    if size(traj{i,1}) < maxNumIter
        plot3(traj{i,1}(:,1),traj{i,1}(:,3),traj{i,1}(:,2),'Color','#EDB120');    % color option: ,'Color','#EDB120'
    end
    hold on;
end
grid on
xlim([-0.397 0.397])
ylim([-0.397 0.397])
zlim([0 1.52])
xlabel('$x$ [m]', 'Interpreter','latex');
ylabel('$z$ [m]', 'Interpreter','latex');
zlabel('$y$ [m]', 'Interpreter','latex');
% view(0,0)
%% Plot trajectories colored with relative velocity for a given injection diameter
figure
% set plot size
set(gcf,'position', [plt_x0,plt_y0,plt_widthRelVel,plt_height])
for i = idx_dTraj:idx_dTraj+24
    if size(traj{i,1}) < maxNumIter
        zRelVelD = sqrt( ( traj{i,1}(:,12) -  traj{i,1}(:,15) ).^2 + ...
                                     ( traj{i,1}(:,13) -  traj{i,1}(:,16) ).^2 + ...
                                     ( traj{i,1}(:,14) -  traj{i,1}(:,17) ).^2  );
        scatter3(traj{i,1}(:,1),traj{i,1}(:,3),traj{i,1}(:,2), 1, zRelVelD)
    end
    hold on;
end
grid on
xlim([-0.397 0.397])
ylim([-0.397 0.397])
zlim([0 1.52])
xlabel('$x$ [m]', 'Interpreter','latex');
ylabel('$z$ [m]', 'Interpreter','latex');
zlabel('$y$ [m]', 'Interpreter','latex');
zCBar = colorbar;
zCBar.Label.String = 'Partocle relative velocity to gas [m/s]';
zCBar.Label.Interpreter = 'latex';
%% Plot trajectories colored with relative velocity for one trajectory in all injection diameters
figure
% set plot size
set(gcf,'position', [plt_x0,plt_y0,plt_widthRelVel,plt_height])
for i = 1:numStream:size(traj,1)
    if size(traj{i,1}) < maxNumIter
        zRelVelD = sqrt( ( traj{i,1}(:,12) -  traj{i,1}(:,15) ).^2 + ...
                                     ( traj{i,1}(:,13) -  traj{i,1}(:,16) ).^2 + ...
                                     ( traj{i,1}(:,14) -  traj{i,1}(:,17) ).^2  );
        scatter3(traj{i,1}(:,1),traj{i,1}(:,3),traj{i,1}(:,2), 1, zRelVelD)
    end
    hold on;
end
grid on
xlim([-0.397 0.397])
ylim([-0.397 0.397])
zlim([0 1.52])
xlabel('$x$ [m]', 'Interpreter','latex');
ylabel('$z$ [m]', 'Interpreter','latex');
zlabel('$y$ [m]', 'Interpreter','latex');
zCBar = colorbar;
zCBar.Label.String = 'Relative velocity magnitude [m/s]';
zCBar.Label.Interpreter = 'latex';
%% Plot statistics of imcomplete trajectory #Restitution 0.99 and 0.7
zz_incomplete = [24 0 1; 22 0 3; 20 1 4; 16 3 6; 18 3 4; 14 2 9; 15 1 9; 12 0 13; 11 1 13];
zz_d = categorical({'22', '27', '32', '37', '42', '47', '52', '57', '60'});
zz_d = reordercats(zz_d,{'22', '27', '32', '37', '42', '47', '52', '57', '60'});
figure
set(gcf,'renderer','Painters')
zz_b = bar(zz_d, zz_incomplete,'stacked');
set(zz_b,'FaceColor','flat');
zz_b(1).CData = [0.4660 0.6740 0.1880]; % complete trajectories: green
zz_b(3).CData = [0.8500 0.3250 0.0980];  % incomplete trajectories: red
zz_b(2).CData = [0 0.4470 0.7410]; % accepted trajectories: blue
ylim([0 35]);
xlabel('Injection diameter [$\mu$m]', 'Interpreter','latex');
ylabel('Number of trajectories [$-$]', 'Interpreter','latex');
legend('Complete','Accepted complete','Incomplete',...
    'Interpreter','latex','Location','best');
%% Plot statistics of imcomplete trajectory #Restitution 0.9 and 0.5
zz_incomplete = [24 1 0; 24 1 0; 21 4 0; 24 1 0; 23 2 0; 25 0 0; 23 2 0; 22 3 0; 16 9 0; 19 6 0];
zz_d = categorical({'19', '22', '27', '32', '37', '42', '47', '52', '57', '60'});
zz_d = reordercats(zz_d,{'19', '22', '27', '32', '37', '42', '47', '52', '57', '60'});
figure
set(gcf,'renderer','Painters')
zz_b = bar(zz_d, zz_incomplete,'stacked');
set(zz_b,'FaceColor','flat');
zz_b(1).CData = [0.4660 0.6740 0.1880]; % complete trajectories: green
zz_b(2).CData = [0.8500 0.3250 0.0980];  % incomplete trajectories: red
zz_b(3).CData = [0 0.4470 0.7410]; % accepted trajectories: blue
ylim([0 35]);
xlabel('Injection diameter [$\mu$m]', 'Interpreter','latex');
ylabel('Number of trajectories [$-$]', 'Interpreter','latex');
legend('Complete','Accepted complete','Incomplete',...
    'Interpreter','latex','Location','best');