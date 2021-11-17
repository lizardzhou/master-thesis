% Use after running dataDistribution
% Need information of outDataPSDg (fine product) and
%                                    outDataPSDpAll (coarse product incl. incomplete trajectories near bottom)
% information includes:
%       (:,1) y-coordinate
%       (:,2) particle diameter at start point
%       (:,3) particle diameter at end point
%       (:,4) trajectory mass flow at start point
%       (:,5) single particle mass at start point
%       (:,6) single particle mass at end point
%       (:,7) trajectory mass flow at end point

% by Xiye Zhou, Oct. 2020
%% mass flow of products and the corresponding end diameter
mFine = [outDataPSDg(:,3), outDataPSDg(:,7)];
mCoarseAll = [outDataPSDpAll(:,3), outDataPSDpAll(:,7)]; % coarse flow incl. accepted
mCoarse = [outDataPSDp(:,3), outDataPSDp(:,7)]; % coarse flow without accepted
%% generate array for calculation
mFlowF = sum(outDataPSDg(:,7)); % total fine flow
mFlowCAll = sum(outDataPSDpAll(:,7)); % total coarse flow incl. accepted
mFlowC = sum(outDataPSDp(:,7)); % total coarse flow without accepted
mFlowAll = mFlowF + mFlowCAll; % total flow incl. accepted
mFlow = mFlowF + mFlowC; % total flow without accepted
fAll = mFlowF / mFlowAll; % fraction of fine product incl. accepted
cAll = 1 - fAll; % fraction of coarse product incl. accepted
f = mFlowF / mFlow; % fraction of fine product without accepted
c = 1 - f; % fraction of coarse product without accepted
dClassG = [0,1:2:15,19:4:35]'*1e-6;
dClassP = [0,20:3:38]'*1e-6;
dMax = 41e-6;
dMin = union(dClassG, dClassP);
dMax = [dMin(2:end); dMax];
% generate array for PSD of separation
CalcArrayAll = zeros(size(dMin)); % incl. accepted
CalcArray = zeros(size(dMin)); % without accepted
CalcArrayAll(:,1) = dMin; % incl. accepted
CalcArrayAll(:,2) = dMax;
CalcArray(:,1) = dMin; % without accepted
CalcArray(:,2) = dMax;
CalcArrayAll(:,3) = (CalcArrayAll(:,2) + CalcArrayAll(:,1)) / 2; % d_m,i
CalcArrayAll(:,4) = CalcArrayAll(:,2) - CalcArrayAll(:,1); % Delta d_i
CalcArray(:,3) = (CalcArray(:,2) + CalcArray(:,1)) / 2; % d_m,i
CalcArray(:,4) = CalcArray(:,2) - CalcArray(:,1); % Delta d_i
%% mass flow of fine product
for i = 1:size(CalcArrayAll,1)
    idx_fineAll = find( (mFine(:,1)>CalcArrayAll(i,1)) & (mFine(:,1)<CalcArrayAll(i,2)) );
    idx_fine = find( (mFine(:,1)>CalcArray(i,1)) & (mFine(:,1)<CalcArray(i,2)) );
    CalcArrayAll(i,5) = sum( mFine(idx_fineAll,2) );
    CalcArray(i,5) = sum( mFine(idx_fine,2) );
end
% mass flow of coarse product
for i = 1:size(CalcArrayAll,1)
    idx_coarseAll = find( (mCoarseAll(:,1)>CalcArrayAll(i,1)) & (mCoarseAll(:,1)<CalcArrayAll(i,2)) );
    idx_coarse = find( (mCoarse(:,1)>CalcArray(i,1)) & (mCoarse(:,1)<CalcArray(i,2)) );
    CalcArrayAll(i,6) = sum( mCoarseAll(idx_coarseAll,2) );
    CalcArray(i,6) = sum( mCoarse(idx_coarse,2) );
end
% omit empty diameter classes
for i = 1:size(CalcArrayAll,1)
    if  CalcArrayAll(i,5) == 0 && CalcArrayAll(i,6) == 0
        CalcArrayAll(i-1,2) = CalcArrayAll(i,2);
    end
    if  CalcArray(i,5) == 0 && CalcArray(i,6) == 0
        CalcArray(i-1,2) = CalcArray(i,2);
    end
end
idx_notEmptyAll = find((CalcArrayAll(:,5)~=0) | (CalcArrayAll(:,6)~=0));
CalcArrayAll = CalcArrayAll(idx_notEmptyAll,:);
idx_notEmpty = find((CalcArray(:,5)~=0) | (CalcArray(:,6)~=0));
CalcArray = CalcArray(idx_notEmpty,:);
%% calculate q_3 for fine and coarse products
CalcArrayAll(:,7) = CalcArrayAll(:,5) ./ (mFlowF * CalcArrayAll(:,4)); % q_fine
CalcArrayAll(:,8) = CalcArrayAll(:,6) ./ (mFlowCAll * CalcArrayAll(:,4)); % q_coarse
CalcArrayAll(:,9) = cAll * CalcArrayAll(:,8) ./ ( fAll * CalcArrayAll(:,7) + cAll * CalcArrayAll(:,8)); % coarse 
CalcArrayAll(:,10) = fAll * CalcArrayAll(:,7) ./ ( fAll * CalcArrayAll(:,7) + cAll * CalcArrayAll(:,8)); % fine
CalcArrayAll(:,11) = fAll * CalcArrayAll(:,7) + cAll * CalcArrayAll(:,8); % all
CalcArray(:,7) = CalcArray(:,5) ./ (mFlowF * CalcArray(:,4)); % q_fine
CalcArray(:,8) = CalcArray(:,6) ./ (mFlowC * CalcArray(:,4)); % q_coarse
CalcArray(:,9) = c * CalcArray(:,8) ./ ( f * CalcArray(:,7) + c * CalcArray(:,8)); % coarse
CalcArray(:,10) = f * CalcArray(:,7) ./ ( f * CalcArray(:,7) + c * CalcArray(:,8)); % fine
CalcArray(:,11) = f * CalcArray(:,7) + c * CalcArray(:,8); % all
%% plot
% separation curve considering accepted complete trajectories
figure
set(gcf,'renderer','Painters')
plot(CalcArrayAll(:,3), CalcArrayAll(:,9), 'o-', CalcArray(:,3), CalcArray(:,9), 'o-');
xlabel('$d_{m,i}$ [m]','Interpreter','latex');
ylabel('$T$ [$-$]','Interpreter','latex');
legend('Complete and accepted complete trajectories','Complete trajectories',...
    'Interpreter','latex','Location','best');
grid on
% q3 for feed, fine and coarse considering accepted complete trajectories
figure
set(gcf,'renderer','Painters')
plot(CalcArrayAll(:,3), CalcArrayAll(:,11), 'o-', CalcArrayAll(:,3), fAll*CalcArrayAll(:,7), 'o-', CalcArrayAll(:,3), cAll*CalcArrayAll(:,8), 'o-');
xlabel('$d_{m,i}$ [m]','Interpreter','latex');
ylabel('$q_{3}$ [m$^{-1}$]','Interpreter','latex');
legend('Feed','Fine product ','Coarse product',...
    'Interpreter','latex','Location','best');
grid on
% q3 for feed, fine and coarse without accepted complete trajectories
figure
set(gcf,'renderer','Painters')
plot(CalcArray(:,3), CalcArray(:,11), 'o-', CalcArray(:,3), f*CalcArray(:,7), 'o-', CalcArray(:,3), c*CalcArray(:,8), 'o-');
xlabel('$d_{m,i}$ [m]','Interpreter','latex');
ylabel('$T(d_{m,i})$ [$-$]','Interpreter','latex');
legend('Feed','Fine','Coarse',...
    'Interpreter','latex','Location','best');
grid on