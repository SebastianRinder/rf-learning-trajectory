%clear variables;

%%%% Plot data with varying size according to reward
% figure
% minMarkerSize = 4;
% maxMarkerSize = 30;
% hold on;
% for k = 1:length(allVals)
% plot(allSamples(k, 1), allSamples(k, 2), '.r', 'MarkerSize', (maxMarkerSize - minMarkerSize) * (allVals(k) - min(allVals)) / max(allVals) + minMarkerSize);
% end
%%%

close all;
%clc;
% set(0, 'DefaultFigureVisible', 'off')

%% seed control
seed = rng;
seed = seed.State(1)
 seed = 2215224515
%  seed = 2103988534; %hard without rotation
% seed = 1165350801 %good for showing local vs global?
% seed = 3268008736;
%  seed = 2420527175 %good showcase for reward switch. bad quad model at iteration ~46
% seed = 3143574292 %local vs global, very good illustration
%seed = 1892178703 %local vs global optim for different kls
%seed = 2345272402 %disconnected
%seed = 3916729757
%seed = 1563644028 %bad m-projection
%seed = 3088878663 %very good m-projection (nbGauss = 50)
rng(seed);
rootPlot = '+static_optimization_algs/plots/';
if(~exist(rootPlot,'dir'))
    mkdir(rootPlot);
end

%% objective function
nbGauss = 25;
muRange = 10;
minSigma = 2;
func = static_optimization_algs.RandomGaussianMixtureFunction(nbGauss, muRange, minSigma);
%%% plotting objective
hnd = figure(1);
func.plot();
funSignature = ['func' num2str(seed) '_' num2str(nbGauss) '_' num2str(muRange) '_' num2str(minSigma)];
fname = [rootPlot 'objective_' funSignature];
hgexport(hnd, [fname '.eps']); %this works better than saveas and print

%% optimization algo
optimizers = {...
%       static_optimization_algs.BayesOptGPML;
%     static_optimization_algs.BayesOptRandom;
%      static_optimization_algs.BayesOpt;
    static_optimization_algs.BayesOptLib;
    };

all_evals = {};
all_signatures = {};

% init distribution and co.
optimizerInput.fun = func;

%samples settings
optimizerInput.functionBounds = [-10 -10; 10 10];
optimizerInput.nbInitSamplesBO = 10;
optimizerInput.nbEvalsBO = 150;
optimizerInput.maxNbSamplesBO = optimizerInput.nbEvalsBO;

% externa paths
addpath('../policysearchtoolbox/NLopt/');
% addpath('../riad_tests_lb/gpml-matlab-v4.0-2016-10-19/');
addpath('../bayesopt/matlab/');
gpStuffPath = 'GPstuff-4.7/';
static_optimization_algs.GP.addGPStuffPath(gpStuffPath);

%GP settings
optimizerInput.hypOptimAfter = 1;
optimizerInput.gpHyperOption = 'MAP';
% optimizerInput.kernelType = 'linear';
optimizerInput.kernelType = 'sexp';
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
    optimizerInput.videoFile = VideoWriter([rootPlot 'policy_search_' all_signatures{end}]);%, 'Uncompressed AVI');
    optimizerInput.videoFile.FrameRate = 4;
    disp(['starting ' all_signatures{end}]);
    tic
    [all_evals{end+1}] = optimizers{k}.optimizeStruct(optimizerInput, func);
    toc
end

%% performance plotting
fname = [rootPlot 'perfOn_' funSignature '_' all_signatures{1}];
hnd = figure(5);
hold on;
for k = 1:length(optimizers)
    plot(all_evals{k});
end
legHnd = legend(all_signatures{:}, 'Location','southeast');
set(legHnd, 'interpreter', 'none');
set(legHnd, 'FontSize', 7);
hgexport(hnd, [fname '.eps']);
