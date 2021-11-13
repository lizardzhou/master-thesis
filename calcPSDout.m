function [d32,dMode,dMedian,calcPSDmat] = calcPSDout(dLowScale,dMax,data,typeProduct)
% Calculate particle size distribution (PSD) of particles escaped as fine or
% coarse products
% Input:
%   dLowScale: defines lower limit of each diameter class (d_min,i) in meter, must be a row vector
%   dMax: defines upper limit of largest diameter class in meter
%   data: particle information given in matrix form, must contain following information:
%       (:,1) y
%       (:,2) particle diameter at start point
%       (:,3) particle diameter at end point
%       (:,4) trajectory mass flow at start point
%       (:,5) single particle mass at start point
%       (:,6) single particle mass at end point
%       (:,7) trajectory mass flow at end point
%   typeProduct: string 'fine' or 'coarse'
% Return:
%   Sauter diameter (d32), mode diameter (dMode), median diameter (dMedian)
%   matrix for calculating PSD (calcPSDmat) containing following information:
%       (:,1) d_min,i
%       (:,2) d_max,i
%       (:,3) d_m,i
%       (:,4) Delta d_i
%       (:,5) mass flow_i
%       (:,6) q_3,i
%       (:,7) summand of Q_3,i
%       (:,8) Q_3,i
%       (:,9) summand of M_-1,3
%       (:,10) summand of M_-3,3
%       (:,11) q_0,i
%       (:,12) summand of q_0,i
%       (:,13) Q_0,i

% by Xiye Zhou, Oct. 2021
%% calculate PSD of mass flow
dataSort = sortrows(data,3); % sort particle diameters ascendingly
calcPSDmat(:,1) = dLowScale'; % lower limit of diameter classes
calcPSDmat(:,2) = [calcPSDmat(2:end,1);dMax]; % upper limit of diameter classes
calcPSDmat(:,3) = (calcPSDmat(:,2) + calcPSDmat(:,1)) / 2; % d_m,i
calcPSDmat(:,4) = calcPSDmat(:,2) - calcPSDmat(:,1); % Delta_d,i
for i = 1:size(calcPSDmat,1)
    outG_idx = find( (dataSort(:,3)>calcPSDmat(i,1)) & (dataSort(:,3)<calcPSDmat(i,2)) );
    calcPSDmat(i,5) = sum(dataSort(outG_idx,7)); % mass flow
end
flowTotal = sum(data(:,7)); % total mass flow
calcPSDmat(:,6) = calcPSDmat(:,5) ./ (flowTotal * calcPSDmat(:,4)); % q_3,i
calcPSDmat(:,7) = calcPSDmat(:,5) / flowTotal; % summand of Q_3,i
calcPSDmat(:,8) = cumsum(calcPSDmat(:,7)); % Q_3,i
calcPSDmat(:,9) = calcPSDmat(:,3).^(-1) .* calcPSDmat(:,4) .* calcPSDmat(:,6); % summand of M_-1,3
d32 = 1 / sum(calcPSDmat(:,9));
[~,idxMax] = max(calcPSDmat(:,6)); % index of the maximum q_3,i
dMode = calcPSDmat(idxMax,3); % mode diameter
dMedian = interp1(calcPSDmat(:,8),calcPSDmat(:,2),0.5); % median diameter
fprintf('Sauter diameter of %s product is %.3g \x03bcm. \n', typeProduct, d32*1e6);
fprintf('Mode diameter of %s product is %.3g \x03bcm. \n', typeProduct, dMode*1e6);
fprintf('Median diameter of %s product is %.3g \x03bcm. \n \n', typeProduct, dMedian*1e6);
%% calculate PSD of number
% q_0,i = q_3,i * d_m,i^(-3) / M_-3,3
calcPSDmat(:,10) = calcPSDmat(:,3).^(-3) .*calcPSDmat(:,4) .* calcPSDmat(:,6); % summand of M_-3,3
calcPSDmat(:,11) = calcPSDmat(:,6) .* calcPSDmat(:,3).^(-3) / sum(calcPSDmat(:,10)); %q_0,i
calcPSDmat(:,12) = calcPSDmat(:,11) .* calcPSDmat(:,4); % summand of Q_0,i
calcPSDmat(:,13) = cumsum(calcPSDmat(:,12)); % Q_0,i
%% plot PSD of mass flow
figure
plot(calcPSDmat(:,3),calcPSDmat(:,6),'o-');
title("Density distribution of mass flow for " + typeProduct + " product");
xlabel('Particle diameter $d_{m,i}$ [m]','Interpreter','latex');
ylabel('$q_3$ [m$^{-1}$]','Interpreter','latex');
grid on;

figure
plot(calcPSDmat(:,2),calcPSDmat(:,8),'o-');
title("Cumulative distribution of mass flow for " + typeProduct + " product");
xlabel('Particle diameter $d_{max,i}$ [m]','Interpreter','latex');
ylabel('$Q_3$ [$-$]','Interpreter','latex');
grid on;
%% plot PSD of number
figure
plot(calcPSDmat(:,3),calcPSDmat(:,11),'o-');
set(gca,'Xscale','log');
grid on
% title("Density distribution of number for " + typeProduct  + " product");
xlabel('Particle diameter $d_{m,i}$ [m]','Interpreter','latex');
ylabel('$q_0$ [m$^{-1}$]','Interpreter','latex');

figure
plot(calcPSDmat(:,2),calcPSDmat(:,13),'o-');
set(gca,'Xscale','log');
grid on
% title("Cumulative distribution of number for " + typeProduct  + " product");
xlabel('Particle diameter $d_{max,i}$ [m]','Interpreter','latex');
ylabel('$Q_0$ [$-$]','Interpreter','latex');
end

