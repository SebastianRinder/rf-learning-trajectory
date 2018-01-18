close all;
clear variables;
rootPlot = '+static_optimization_algs/plots/';
if(~exist(rootPlot,'dir'))
    mkdir(rootPlot);
end

%% seed initialization
seed = rng;
seed = seed.State(1)
seed = 1299576886
rng(seed);

%% objective function
nbGauss = 25;
muRange = 10;
minSigma = 2;
numSamples = 50;

func = static_optimization_algs.RandomGaussianMixtureFunction(nbGauss, muRange, minSigma);
%%% plotting objective
hnd = figure;
func.plot();
funSignature = ['func' num2str(seed) '_' num2str(nbGauss) '_' num2str(muRange) '_' num2str(minSigma) '_' num2str(numSamples)];
fname = [rootPlot 'objective_' funSignature];

%% sample some points from a Gaussian
seed = 129957
rng(seed);
mu = [-1 -2];
covC = [.25 -.2; -.2 1] * 5;
% mu = [-1 -3.5];
% covC = [.25 -.4; -.4 1] * 5;

dataDis = static_optimization_algs.Normal(mu, covC);
data = dataDis.getSamples(numSamples);
evals = func.eval(data);
hold on;
optPos = [-1.79, 0.47];
opt2Pos = [-.7, -3.9];

plot(optPos(1), optPos(2), '*r', 'MarkerSize', 20);
plot(opt2Pos(1), opt2Pos(2), '*b', 'MarkerSize', 20);
axis([-4,2, -8,4]);

hgexport(hnd, [fname '.eps']);

hnd = figure;
dataDis.plot();
hold on;
plot(optPos(1), optPos(2), '*r', 'MarkerSize', 20);
plot(opt2Pos(1), opt2Pos(2), '*b', 'MarkerSize', 20);
plot(data(:, 1), data(:, 2), '*g');
axis([-4,2, -8,4]);
hgexport(hnd, [fname '_dist.eps']);
%% quad model
% regul = 1;
% for k = 1:5
%     subplot(1,6,k+1);
%     [~, quadModel] = static_optimization_algs.Quadratic.learnQuadModel(data, evals, regul);
%     regul = regul * 10;
%     static_optimization_algs.Quadratic.plotQuad(quadModel, [-1, -2], [3, 6]);
%     hold on;
%     plot(optPos(1), optPos(2), '*r', 'MarkerSize', 20);
%     modelOpt = -quadModel.A \ quadModel.w;
%     plot(modelOpt(1), modelOpt(2), '*m', 'MarkerSize', 20);
%     axis([-4,2, -8,4]);
% end

hnd = figure;
[~, quadModel] = static_optimization_algs.Quadratic.learnQuadModel(data, evals, 1e-6);
static_optimization_algs.Quadratic.plotQuad(quadModel, [-1, -2], [3, 6]);
hold on;
plot(optPos(1), optPos(2), '*r', 'MarkerSize', 20);
plot(opt2Pos(1), opt2Pos(2), '*b', 'MarkerSize', 20);
modelOpt = -quadModel.A \ quadModel.w;
plot(modelOpt(1), modelOpt(2), '*m', 'MarkerSize', 20);
axis([-4,2, -8,4]);
hgexport(hnd, [fname '_quad.eps']);

%metric test
cC = chol(covC);
invCC = cC \ eye(size(cC));

x = data * invCC;
max(abs(pdist(x) - pdist(data, 'mahalanobis', covC)))

%% gp regression
static_optimization_algs.GP.addGPStuffPath('GPstuff-4.7/');

x = data * invCC;
y = evals - mean(evals);

% likelihood and covraiance function
lik = lik_gaussian('sigma2', 0.2^2);
gpcf = gpcf_sexp('lengthScale', 1, 'magnSigma2', 0.2^2);
% hyper-param priors
pn = prior_logunif();
lik = lik_gaussian(lik,'sigma2_prior', pn);
pl = prior_unif();
pm = prior_sqrtunif();
gpcf = gpcf_sexp(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
gp = gp_set('lik', lik, 'cf', gpcf);

% MCMC (slice sampling)
disp('MCMC integration over the parameters')
[gp_rec,g,opt] = gp_mc(gp, x, y, 'nsamples', 420);
gp_rec = thin(gp_rec, 20, 4); %delete first 20 samples and take every second sample

% plotting the GP
hnd = figure;
static_optimization_algs.GP.plotGP(gp_rec, x, y, [-1, -2], [3, 6], [], invCC);
hold on;
plot(optPos(1), optPos(2), '*r', 'MarkerSize', 20);
plot(opt2Pos(1), opt2Pos(2), '*b', 'MarkerSize', 20);
axis([-4,2, -8,4]);
hgexport(hnd, [fname '_gp.eps']);

hnd = figure;
samplesArgmax = zeros(length(gp_rec.e), 2);
for k = 1:length(gp_rec.e)
    gp.lik.sigma2 = gp_rec.lik.sigma2(k);
    gp.cf{1}.lengthScale = gp_rec.cf{1}.lengthScale(k);
    gp.cf{1}.magnSigma2 = gp_rec.cf{1}.magnSigma2(k);
    newDist = static_optimization_algs.Normal(dataDis.mu, 1 * dataDis.getCov);
    samplesArgmax(k, :) = static_optimization_algs.GP.localArgmaxSample(gp, x, y, newDist, 800, [], invCC);
end
plot(samplesArgmax(:, 1), samplesArgmax(:, 2), '*g');
hold on;
amaxApprox = static_optimization_algs.Normal([0 0], eye(2));
amaxApprox = amaxApprox.wmle(samplesArgmax);
plot(optPos(1), optPos(2), '*r', 'MarkerSize', 20);
plot(amaxApprox.mu(1), amaxApprox.mu(2), '*m', 'MarkerSize', 20);
plot(opt2Pos(1), opt2Pos(2), '*b', 'MarkerSize', 20);
amaxApprox.plot
axis([-4,2, -8,4]);
hgexport(hnd, [fname '_argmax.eps']);

figure
idx = randperm(length(gp_rec.cf{1}.magnSigma2));
for k = 1:5
    subplot(1, 5, k);
    gp.lik.sigma2 = gp_rec.lik.sigma2(idx(k));
    gp.cf{1}.lengthScale = gp_rec.cf{1}.lengthScale(idx(k));
    gp.cf{1}.magnSigma2 = gp_rec.cf{1}.magnSigma2(idx(k));
    static_optimization_algs.GP.plotGP(gp, x, y, [-1, -2], [3, 6], [], invCC);
end


