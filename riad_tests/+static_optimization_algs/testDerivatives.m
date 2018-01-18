%test derivatives
clear variables;
close all;

seed = rng;
seed = seed.State(1)
%seed = 4012388177 
%seed = 2345272402 %disconnected
%seed = 443676256
seed = 3671025236
rng(seed);

dim = 1000;

nbGauss = 25;
muRange = 10;
minSigma = 2;
func = static_optimization_algs.RandomGaussianMixtureFunction(nbGauss, muRange, minSigma, dim);

% QuadFunction Test
% rangeR = 10;
% r = rangeR * (ones(dim, 1) - 2 * rand(dim, 1));
% a0 = zeros(dim * dim, 1);
% [y, dx] = static_optimization_algs.Test.do(a0, r)
% diff = static_optimization_algs.GradientCheck.do(@(a) static_optimization_algs.Test.do(a, r), a0, 1e-4)

rangeMu = 0;
minVar = 1;
covC = eye(dim) * minVar;
covC = ones(dim, dim) - 2 / sqrt(dim) * rand(dim, dim);
covC = covC + diag(minVar * ones(dim, 1));
covC = covC' * covC;

eyeDim = eye(dim);
% disp('inv(covC)');
% tic
% invCovC = inv(covC);
% toc
% 
disp('covC\eyeDim');
tic
invCovCDiv = covC\eyeDim;
toc

disp('chol(covC)');
tic
cholCovC = chol(covC);
toc

% disp('detCovC');
% tic
% detCovC = det(covC);
% toc
% 
% disp('detcholc');
% tic
% detCholC = det(cholCovC);
% toc
% disp('detcholcByDiag')
% tic
% detCholCDiag = prod(diag(cholCovC));
% toc
% disp('inv(cholc)');
% tic
% invCholC = inv(cholCovC);
% toc
disp('inv(cholc)');
tic
invCholC = cholCovC\eyeDim;
toc
% disp('inv cholc sparse')
% sCholC = sparse(cholCovC);
% tic
% sInvCholC = sCholC\eyeDim;
% toc
% 
% disp('cholc * covC');
% tic
% prod2 = cholCovC * covC;
% toc
% disp('scholc * covC')
% tic
% prod = sCholC * covC;
% toc


% disp('eig(cholc)');
% tic
% [V,D,W] = eig(covC);
% toc


% 
% disp('invcholc * invcholc''');
% tic
% precision = invCholC * invCholC';
% toc
% 

mu = rangeMu * (ones(dim, 1) - 2 * rand(dim, 1));
normal = static_optimization_algs.Normal(mu, covC);

nbPoints = 100;
precision = inv(covC);
params0 = [mu; precision(tril(ones(dim)) == 1)];
thetas = normal.getSamples(nbPoints);
probas = exp(normal.getLogProbas(thetas));
rewards = func.eval(thetas);
%rewards = rewards / max(rewards);
% optimization
eta = .01;
% lbVar = 1e-2;
% ubVar = 1/lbVar;
% lbParam = -inf * ones(size(params0));
% ubParam =  inf * ones(size(params0));
% varIndexes = eye(dim);
% varIndexes = [zeros(dim, 1); varIndexes(tril(ones(dim)) == 1)];
% lbParam(varIndexes == 1) = lbVar * ones(dim, 1);
% ubParam(varIndexes == 1) = ubVar * ones(dim, 1);
options =  optimset('GradObj','on', 'DerivativeCheck', 'on', 'Display', 'on');
% options.TolFun = 0.1;

%[diff, numGrad] = static_optimization_algs.GradientCheck.do(@(params) static_optimization_algs.RepsBanditsGradReverse.Lagrangian(params, covC, det(covC), mu', thetas, rewards, probas, 0.1, eta), params0, 1e-4)

% [yBefore, dxBefore] = static_optimization_algs.IProjection.eval(params0, dim, precision, mu', thetas, rewards, probas, eta);
% tic
% params = fmincon(@(params) static_optimization_algs.IProjection.eval(params, dim, precision, mu', thetas, rewards, probas, eta), params0, [], [], [], [], [], [], [], options);
% toc
% [yAfter, dxAfter] = static_optimization_algs.IProjection.eval(params, dim, precision, mu', thetas, rewards, probas, eta);
% [diff, numGrad] = static_optimization_algs.GradientCheck.do(@(paramsv) static_optimization_algs.IProjection.eval(paramsv, dim, precision, mu', thetas, rewards, probas, eta), params, 1e-4)
% 
% params0 = params;
% tic
% params = fmincon(@(params) static_optimization_algs.IProjection.eval(params, dim, precision, mu', thetas, rewards, probas, eta), params0, [], [], [], [], [], [], [], options);
% toc
% [yAfter, dxAfter] = static_optimization_algs.IProjection.eval(params, dim, precision, mu', thetas, rewards, probas, eta);

lbParams = -inf * ones(size(params0));
lbParams(dim+1:end) = 0;

[yBefore, dxBefore] = static_optimization_algs.IProjection.eval_chol(params0, dim, precision, mu', thetas, rewards, probas, eta);
tic
params = fmincon(@(params) static_optimization_algs.IProjection.eval_chol(params, dim, precision, mu', thetas, rewards, probas, eta), params0, [], [], [], [], lbParams, [], [], options);
toc
[yAfter, dxAfter] = static_optimization_algs.IProjection.eval_chol(params, dim, precision, mu', thetas, rewards, probas, eta);

params0 = params;
tic
params = fmincon(@(params) static_optimization_algs.IProjection.eval_chol(params, dim, precision, mu', thetas, rewards, probas, eta), params0, [], [], [], [], lbParams, [], [], options);
toc
[yAfter, dxAfter] = static_optimization_algs.IProjection.eval_chol(params, dim, precision, mu', thetas, rewards, probas, eta);


%[diff, numGrad] = static_optimization_algs.GradientCheck.do(@(paramsv) static_optimization_algs.IProjection.eval(paramsv, dim, precision, mu', thetas, rewards, probas, eta), params, 1e-4)
% tic 
% [yBefore, dxBefore] = static_optimization_algs.RepsBanditsGradReverse.Lagrangian(params0, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta);
% toc
% 
% [diff, numGrad] = static_optimization_algs.GradientCheck.do(@(paramsv) static_optimization_algs.RepsBanditsGradReverse.Lagrangian(paramsv, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta), params0, 1e-4)
% tic
% params = fmincon(@(params) static_optimization_algs.RepsBanditsGradReverse.Lagrangian(params, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta), params0, [], [], [], [], [], [], [], options);
% toc
% [yAfter, dxAfter] = static_optimization_algs.RepsBanditsGradReverse.Lagrangian(params, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta);
% 
% params0 = params;
% [diff, numGrad] = static_optimization_algs.GradientCheck.do(@(paramsv) static_optimization_algs.RepsBanditsGradReverse.Lagrangian(paramsv, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta), params0, 1e-4)
% 
% tic
% params = fmincon(@(params) static_optimization_algs.RepsBanditsGradReverse.Lagrangian(params, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta), params0, [], [], [], [], [], [], [], options);
% toc
% [yAfter, dxAfter] = static_optimization_algs.RepsBanditsGradReverse.Lagrangian(params, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta);


% [diff, numGrad] = static_optimization_algs.GradientCheck.do(@(paramsv) static_optimization_algs.RepsBanditsGradReverse.Lagrangian(paramsv, diag(diag(covC)), prod(diag(covC)), mu', thetas, rewards, probas, 0.1, eta), params, 1e-4)

% plotting
figure('Position', [200, 0, 1600, 900]);

subplot(1,3,1);
hold on;
func.plot();

rewardPlot = rewards ./ eta;
minMarkerSize = 1;
maxMarkerSize = 20;
aPlot = (maxMarkerSize - minMarkerSize) / (max(rewardPlot) - min(rewardPlot));
bPlot = minMarkerSize - aPlot * min(rewardPlot);
for k = 1:nbPoints
    plot(thetas(k, 1), thetas(k, 2), '.r', ...
        'MarkerSize', floor(aPlot * rewardPlot(k) + bPlot));
end
normal.plot();

subplot(1,3,2);
hold on;
func.plot();
plot(thetas(:, 1), thetas(:, 2), '.b');
normalMle = normal.wmle(thetas, exp((rewards-max(rewards))/eta));
normalMle.plot();

subplot(1,3,3);
hold on;
func.plot();
plot(thetas(:, 1), thetas(:, 2), '.b');
muOptim = params(1:dim);
covOptim = zeros(dim);
covOptim(tril(ones(dim)) == 1) = params(dim+1:end);
covOptim = covOptim + covOptim' - diag(diag(covOptim));
covOptim = inv(covOptim);
normalOptim = static_optimization_algs.Normal(muOptim, covOptim);
normalOptim.plot();
[y, dx] = static_optimization_algs.IProjection.eval(params, dim, precision, mu', thetas, rewards, probas, eta);



% [diff, numGrad] = static_optimization_algs.GradientCheck.do(@(params) static_optimization_algs.Test.WeightedGauss(params, dim, thetas, weights), params0, 1e-6);
% diff 
% dcovNum = zeros(dim);
% dcovNum(tril(ones(dim)) == 1) = numGrad(dim+1:end);
% 
% 
