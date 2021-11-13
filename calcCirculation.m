%% relative velocity magnitude
figure
for idxTraj1 = 176:200
    relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                             ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                             ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
    tau1 = traj{idxTraj1,1}(:,4);
    plot(tau1, relVel_Stk,'Color','#0072BD');
    hold on
end
set(gca,'Xscale','log');
xlim([1e-6 1e3]);
plt1 = plot(tau1, relVel_Stk,'Color','#0072BD');
hold on
for idxTraj1 = 276:300
    relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                             ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                             ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
    tau1 = traj{idxTraj1,1}(:,4);
    plot(tau1, relVel_Stk,'Color','#D95319');
    hold on
end
plt2 = plot(tau1, relVel_Stk,'Color','#D95319');
hold on
for idxTraj1 = 376:400
    relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                             ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                             ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
    tau1 = traj{idxTraj1,1}(:,4);
    plot(tau1, relVel_Stk,'Color','#7E2F8E');
    hold on
end
plt3 = plot(tau1, relVel_Stk,'Color','#7E2F8E');
legend([plt1,plt2,plt3],'22','42','60')
%% Relative Reynolds number
viscosityG = 1.72e-5;
densityG = 1.039;
figure
for idxTraj1 = 176:200
    d =  traj{idxTraj1,1}(:,6);
    relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                             ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                             ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
    ReRel = d .* relVel_Stk * densityG / viscosityG; % relative Reynolds number
    tau1 = traj{idxTraj1,1}(:,4);
    plot(tau1, ReRel,'Color','#0072BD');
    hold on
end
set(gca,'Xscale','log');
plt1 = plot(tau1, ReRel,'Color','#0072BD');
hold on
for idxTraj1 = 276:300
    d =  traj{idxTraj1,1}(:,6);
    relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                             ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                             ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
    ReRel = d .* relVel_Stk * densityG / viscosityG;
    tau1 = traj{idxTraj1,1}(:,4);
    plot(tau1, ReRel,'Color','#D95319');
    hold on
end
plt2 = plot(tau1, ReRel,'Color','#D95319');
hold on
for idxTraj1 = 376:400
    d =  traj{idxTraj1,1}(:,6);
    relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                             ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                             ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
    ReRel = d .* relVel_Stk * densityG / viscosityG;
    tau1 = traj{idxTraj1,1}(:,4);
    plot(tau1, ReRel,'Color','#7E2F8E');
    hold on
end
plt3 = plot(tau1, ReRel,'Color','#7E2F8E');
legend([plt1,plt2,plt3], '22', '42', '60')
%% Stokes number
maxNumIter = 150e3;
plt_x0 = 100;
plt_y0 = 100;
plt_width = 900;
plt_height = 350;
figure
set(gcf,'position', [plt_x0,plt_y0,plt_width,plt_height])
set(gcf,'renderer','Painters')
for idxTraj1 = 176:200 % 22 microns
    if size(traj{idxTraj1,1},1) < maxNumIter
        d =  traj{idxTraj1,1}(:,6);
        velG = sqrt(  traj{idxTraj1,1}(:,15).^2 + traj{idxTraj1,1}(:,16).^2 + traj{idxTraj1,1}(:,17).^2 );
        relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                                       ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                                       ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
        ReRel = d .* relVel_Stk * densityG / viscosityG;
        dragCoeff = zeros( size(traj{idxTraj1,1},1), 1);
        for jj = 1:size(traj{idxTraj1,1},1) 
               if ReRel(jj) <0.1
                    dragCoeff(jj) = 24 / ReRel(jj);
               else
                    if ReRel(jj) < 1
                        dragCoeff(jj) = 22.73 ./ ReRel(jj) + 0.0903 ./ ReRel(jj).^2 + 3.69;
                    else
                        if ReRel(jj) < 10
                            dragCoeff(jj) = 29.1667 ./ ReRel(jj) - 3.8889 ./ ReRel(jj).^2 + 1.222;
                        else
                            if ReRel(jj) < 100
                                dragCoeff(jj) = 46.5 ./ ReRel(jj) - 116.67 ./ ReRel(jj).^2 + 0.6167;
                            else
                                dragCoeff(jj) = 98.33 ./ ReRel(jj) - 2778 ./ ReRel(jj).^2 + 0.3644;
                            end
                        end
                    end
               end
        end
        relaxTau = densityG .* d.^2 * 24 ./ (18 * viscosityG .* dragCoeff .* ReRel);
        Stk = relaxTau .* velG ./ d;
        tau1 = traj{idxTraj1,1}(:,4);
        plot(tau1, Stk,'Color','#0072BD');
        hold on
    end
end
set(gca,'Xscale','log');
xlim([1e-6 1e3]);
ylim([0 2])
yline(1,'--');
plt1 = plot(tau1, Stk,'Color','#0072BD');
grid on
hold on
xlabel('Residence time [s]','Interpreter','latex');
ylabel('Stokes number [$-$]','Interpreter','latex');
legend(plt1,'Injection diameter 22\,$\mu$m', 'Interpreter','latex')

figure
set(gcf,'position', [plt_x0,plt_y0,plt_width,plt_height])
set(gcf,'renderer','Painters')
idxc = 1;
for idxTraj1 =  276:300 % 42 microns
    if size(traj{idxTraj1,1},1) < maxNumIter
        d =  traj{idxTraj1,1}(:,6);
        velG = sqrt(  traj{idxTraj1,1}(:,15).^2 + traj{idxTraj1,1}(:,16).^2 + traj{idxTraj1,1}(:,17).^2 );
        relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                                       ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                                       ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
        ReRel = d .* relVel_Stk * densityG / viscosityG;
        dragCoeff = zeros( size(traj{idxTraj1,1},1), 1);
        for jj = 1:size(traj{idxTraj1,1},1) 
               if ReRel(jj) <0.1
                    dragCoeff(jj) = 24 / ReRel(jj);
               else
                    if ReRel(jj) < 1
                        dragCoeff(jj) = 22.73 ./ ReRel(jj) + 0.0903 ./ ReRel(jj).^2 + 3.69;
                    else
                        if ReRel(jj) < 10
                            dragCoeff(jj) = 29.1667 ./ ReRel(jj) - 3.8889 ./ ReRel(jj).^2 + 1.222;
                        else
                            if ReRel(jj) < 100
                                dragCoeff(jj) = 46.5 ./ ReRel(jj) - 116.67 ./ ReRel(jj).^2 + 0.6167;
                            else
                                dragCoeff(jj) = 98.33 ./ ReRel(jj) - 2778 ./ ReRel(jj).^2 + 0.3644;
                            end
                        end
                    end
               end
        end
        relaxTau = densityG * d.^2 * 24 ./ (18 * viscosityG .* dragCoeff .* ReRel);
        Stk = relaxTau .* velG ./ d;
        tau1 = traj{idxTraj1,1}(:,4);
        plot(tau1, Stk, 'Color','#D95319'); %,'Color','#D95319'    %, 'Color', stk_cmap(idxc,:)
        hold on
    end
end
fprintf('\n');
plt2 = plot(tau1, Stk,'Color','#D95319');
set(gca,'Xscale','log');
xlim([1e-6 1e3]);
ylim([0 2]);
yline(1,'--');
grid on
xlabel('Residence time [s]','Interpreter','latex');
ylabel('Stokes number [$-$]','Interpreter','latex')
legend(plt2,'Injection diameter 42\,$\mu$m', 'Interpreter','latex')

figure
set(gcf,'position', [plt_x0,plt_y0,plt_width,plt_height])
set(gcf,'renderer','Painters')
for idxTraj1 = 376:400 % 400 microns
    if size(traj{idxTraj1,1},1) < maxNumIter
        d =  traj{idxTraj1,1}(:,6);
        velG = sqrt(  traj{idxTraj1,1}(:,15).^2 + traj{idxTraj1,1}(:,16).^2 + traj{idxTraj1,1}(:,17).^2 );

        relVel_Stk = sqrt( ( traj{idxTraj1,1}(:,12) -  traj{idxTraj1,1}(:,15) ).^2 + ...
                                 ( traj{idxTraj1,1}(:,13) -  traj{idxTraj1,1}(:,16) ).^2 + ...
                                 ( traj{idxTraj1,1}(:,14) -  traj{idxTraj1,1}(:,17) ).^2  );
        ReRel = d .* relVel_Stk * densityG / viscosityG;
        dragCoeff = zeros( size(traj{idxTraj1,1},1), 1);
        for jj = 1:size(traj{idxTraj1,1},1) 
               if ReRel(jj) <0.1
                    dragCoeff(jj) = 24 / ReRel(jj);
               else
                    if ReRel(jj) < 1
                        dragCoeff(jj) = 22.73 ./ ReRel(jj) + 0.0903 ./ ReRel(jj).^2 + 3.69;
                    else
                        if ReRel(jj) < 10
                            dragCoeff(jj) = 29.1667 ./ ReRel(jj) - 3.8889 ./ ReRel(jj).^2 + 1.222;
                        else
                            if ReRel(jj) < 100
                                dragCoeff(jj) = 46.5 ./ ReRel(jj) - 116.67 ./ ReRel(jj).^2 + 0.6167;
                            else
                                dragCoeff(jj) = 98.33 ./ ReRel(jj) - 2778 ./ ReRel(jj).^2 + 0.3644;
                            end
                        end
                    end
               end
        end
        relaxTau = densityG .* d.^2 * 24 ./ (18 * viscosityG .* dragCoeff .* ReRel);
        Stk = relaxTau .* velG ./ d;
        tau1 = traj{idxTraj1,1}(:,4);
        plot(tau1, Stk,'Color','#7E2F8E');
        hold on
    end
end
fprintf('\n');
set(gca,'Xscale','log');
xlim([1e-6 1e3]);
ylim([0 2])
yline(1,'--');
% plt1 = plot(tau1, Stk,'Color','#0072BD');
grid on
hold on
xlabel('Residence time [s]','Interpreter','latex');
ylabel('Stokes number [$-$]','Interpreter','latex')
plt3 = plot(tau1, Stk,'Color','#7E2F8E'); %,'LineWidth',5
legend(plt3,'Injection diameter 60\,$\mu$m', 'Interpreter','latex')
% legend([plt1,plt2,plt3],'22','42','60')