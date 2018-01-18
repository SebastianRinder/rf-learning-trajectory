clear variables;

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
rng('default');
rng('shuffle');
seed = rng;
seed = seed.State(1)
% seed = 837594802
% seed = 59823194
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
% funsAll = benchmarks('handles');
% k = 9;
% func = static_optimization_algs.CocoWrapper(funsAll{k}, ['f' num2str(k)], true)
%%% plotting objective
hnd = figure(1);
func.plot();
funSignature = ['func' num2str(seed) '_' num2str(nbGauss) '_' num2str(muRange) '_' num2str(minSigma)];
fname = [rootPlot 'objective_' funSignature];
hgexport(hnd, [fname '.eps']); %this works better than saveas and print

%% optimization algo
optimizers = {...
      %static_optimization_algs.RepsBandits; ...
%     static_optimization_algs.RepsBanditsIProjection; ...
%     static_optimization_algs.RepsBanditsIProjectionChol; ...
    %     static_optimization_algs.RepsBanditsDecay, ...
%     static_optimization_algs.RepsTrueGrid, ...
    % static_optimization_algs.RepsTrueGaussian,  ...
    %static_optimization_algs.RepsTrueGaussianIProjection; ...
%       static_optimization_algs.RepsBanditsGrad; ...
%       static_optimization_algs.RepsBanditsGradNonParamIS; ...
      %      static_optimization_algs.RepsBanditsGradHoeffding; ...
%     static_optimization_algs.AdaptiveMore;

%     static_optimization_algs.AdaptiveMoreMW_Rotate;
%     static_optimization_algs.AdaptiveMoreWithModelSwitch;
%   
%     static_optimization_algs.AdaptiveMoreAndersonWithModelSwitch;
%     static_optimization_algs.AdaptiveMoreStudent_Rotate;
%    static_optimization_algs.AdaptiveMoreStudentEntropy;
%     static_optimization_algs.AdaptiveMoreStudent;
%     static_optimization_algs.LocalBayes;
    static_optimization_algs.DensityWeightedBO;
%     static_optimization_algs.CMAESWrapper;
%     static_optimization_algs.More;

%     static_optimization_algs.ConservativeAdaptiveMore;
%      static_optimization_algs.AdaptiveMoreIProj;
      %     static_optimization_algs.RepsBanditsGradChol; ...
%     static_optimization_algs.RepsBanditsGradReverse; ...
    };
%!!! change logDetP and how KL is computed using it

all_perfs = {};
all_kls = {};
all_signatures = {};

% init distribution and co.
optimizerInput.fun = func;
optimizerInput.initVar = 10;
mu = [0 0]; covC = eye(2) * optimizerInput.initVar;
optimizerInput.initDistrib = static_optimization_algs.Normal(mu, covC);
optimizerInput.minEpsiKL = .01;
optimizerInput.epsiKL = .05;
optimizerInput.entropyReduction = .05;
optimizerInput.nbSamplesPerIter = 10;
optimizerInput.nbInitSamples = 5;
optimizerInput.nbIter = 30;
optimizerInput.maxEvals = 500; %cmaes wrapper only depends on this.
%optimizerInput.maxIterReuse = optimizerInput.nbIter;
optimizerInput.maxIterReuse = 30;
optimizerInput.sampleDecay = .8;
optimizerInput.nbBinPerDim = 100;

optimizerInput.useImportanceSampling = 1;

optimizerInput.confidenceMW = .1;
optimizerInput.confidenceStudent = .1;

optimizerInput.confidenceCMW = .5;
optimizerInput.regularization = 1e-6;

%local bayes settings
optimizerInput.gpStuffPath = 'GPstuff-4.7/';
optimizerInput.nbPStarSamples = 100;
optimizerInput.minProbaReuse = inf;
% optimizerInput.gpHyperOption = 'MCMC';
optimizerInput.gpHyperOption = 'MAP';

optimizerInput.samplingOption = 'Acquisition';
%optimizerInput.samplingOption = 'Random';

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
    optimizerInput.videoFile = VideoWriter([rootPlot 'policy_search_' all_signatures{end}]);%, 'Uncompressed AVI');
    optimizerInput.videoFile.FrameRate = 4;
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

