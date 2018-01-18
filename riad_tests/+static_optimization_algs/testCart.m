close all;
%% seed control
rng('default');
rng('shuffle');
seed = rng;
% seed = seed.State(1)
% seed = 596266919
rng(seed);
rootPlot = '+static_optimization_algs/plots/';
if(~exist(rootPlot,'dir'))
    mkdir(rootPlot);
end

%% objective function
addpath('cartPole');
func = cartPole();
%%% plotting objective
funSignature = ['cartpole'];
%% optimization algo
optimizers = {...
    static_optimization_algs.DensityWeightedBO;
    };

all_perfs = {};
all_kls = {};
all_signatures = {};

% init distribution and co.
optimizerInput.fun = func;
optimizerInput.initVar = 1.5;
mu = zeros(1, func.opts.dim);
covC = eye(func.opts.dim) * optimizerInput.initVar;
optimizerInput.initDistrib = static_optimization_algs.Normal(mu, covC);
optimizerInput.epsiKL = .05;
optimizerInput.entropyReduction = .05;
optimizerInput.nbSamplesPerIter = 6;
optimizerInput.nbIter = 100;
optimizerInput.maxIterReuse = 30;

%local bayes settings
optimizerInput.gpStuffPath = 'GPstuff-4.7/';
optimizerInput.minProbaReuse = inf;
% optimizerInput.gpHyperOption = 'MCMC';
optimizerInput.gpHyperOption = 'MAP';

optimizerInput.samplingOption = 'Acquisition';

% optimizerInput.kernelType = 'linear';
% optimizerInput.kernelType = 'sexp';
optimizerInput.kernelType = 'matern52';

% optimizerInput.featureFunction = @(x) [ones(size(x, 1), 1) x];
% optimizerInput.featureFunction = @(data) static_optimization_algs.Quadratic.getQuadFeatures(data, 1, 1);
optimizerInput.featureFunction = [];

optimizerInput.featureName = 'noFeature';
% optimizerInput.featureName = 'quadratic';

% optimizerInput.yCenteringType = 'min';
optimizerInput.yCenteringType = 'mean';
% optimizerInput.yCenteringType = 'max';

optimizerInput.videoFile = [];

seed = rng;
seedStartOpt = seed.State(2)

for k = 1:length(optimizers)
    rng(seedStartOpt);
    all_signatures{end+1} = [optimizers{k}.getSignature(optimizerInput) '_' funSignature];
    disp(['starting ' all_signatures{end}]);
    tic
    [all_perfs{end+1}] = optimizers{k}.optimizeStruct(optimizerInput, func);
    toc
end

%% performance plotting
fname = [rootPlot 'perfOn_' funSignature '_' all_signatures{1}];
hnd = figure(5);
hold on;
for k = 1:length(optimizers)
    plot(all_perfs{k});
    perfs = all_perfs{k};
end
legHnd = legend(all_signatures{:}, 'Location','southeast');
set(legHnd, 'interpreter', 'none');
set(legHnd, 'FontSize', 7);
hgexport(hnd, [fname '.eps']);

%% kl plotting
% fname = [rootPlot 'klOn_' funSignature '_' all_signatures{1}];
% hnd = figure(6);
% for k = 1:length(optimizers)
%     semilogy(all_perfs{k}(:, 1), all_kls{k});
%     hold on;
% end
% hold off;
% legHnd = legend(all_signatures{:}, 'Location','southeast');
% set(legHnd, 'interpreter', 'none');
% set(legHnd, 'FontSize', 7);
% hgexport(hnd, [fname '.eps']); 

