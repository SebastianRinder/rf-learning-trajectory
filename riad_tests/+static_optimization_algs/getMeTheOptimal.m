close all;
clear variables;
load decreasedReward_solvedChol_solvedOptim;

% searching for better eta/gamma
% obj.updateModel([]);
% obj.updateModel([]);
% obj.updateModel([]);
% 
% nbReOptim = 5;
% for k = 1:nbReOptim
%     nbPOpt = 50;
%     limitEta = [obj.eta - obj.eta / 1.5; obj.eta + obj.eta / 1.5];
%     limitGamma = [obj.gammA - obj.gammA / 1.5; obj.gammA + obj.gammA / 1.5];
% 
%     lineEta = limitEta(1):(limitEta(2)-limitEta(1))/nbPOpt:limitEta(2);
%     lineGamma = limitGamma(1):(limitGamma(2)-limitGamma(1))/nbPOpt:limitGamma(2);
% 
%     [xEta, xGamma] = meshgrid(lineEta, lineGamma);
%     y = zeros(size(xEta(:)));
%     for i = 1:length(y)
%         obj.eta = xEta(i);
%         obj.gammA = xGamma(i);
%         try
%             obj.updatePolicy(policy);
%             y(i) = obj.rewardAfterUpdate;
%         catch
%             y(i) = -inf;
%         end
%     end
%     figure;
%     h = surf(lineEta, lineGamma, reshape(y, length(lineGamma), length(lineEta)));
%     set(h,'LineStyle','none');
%     [val, argM] = max(y);
%     obj.eta = xEta(argM);
%     obj.gammA = xGamma(argM);
% end

% policy for mean state
policyMeanState = static_optimization_algs.Normal(obj.bias + obj.K * obj.muState,...
    obj.Q + obj.K * obj.SigmaState * obj.K');

% plotting the reward for the average state
nbPoint = 1000;
[limitInf, limitSup] = policyMeanState.getLimits;
areaLimit = 2;
limitInf = limitInf * areaLimit;
limitSup = limitSup * areaLimit;
line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
[x1, x2] = meshgrid(line1, line2);
x = [x1(:) x2(:)];
y = dot(x * obj.Raa, x, 2) + x * obj.Ras * obj.muState + x * obj.ra;
h = surf(line1, line2, reshape(y, length(line2), length(line1)));
set(h,'LineStyle','none');
hold on;
policyMeanState.plot();
[i, v] = max(y);
plot(x(v, 1), x(v, 2), 'xr');
muA = policyMeanState.getMu;
plot(muA(1), muA(2), 'xb');

% update policy and replot
obj.updatePolicy(policy);
figure; 
hold on;
obj.Q = policy.getCovariance();
obj.bias = policy.bias;
obj.K = policy.weights;
policyMeanState = static_optimization_algs.Normal(obj.bias + obj.K * obj.muState,...
    obj.Q + obj.K * obj.SigmaState * obj.K');
h = surf(line1, line2, reshape(y, length(line2), length(line1)));
set(h,'LineStyle','none');
hold on;
policyMeanState.plot();
[i, v] = max(y);
plot(x(v, 1), x(v, 2), 'xr');
muA = policyMeanState.getMu;
plot(muA(1), muA(2), 'xb');




