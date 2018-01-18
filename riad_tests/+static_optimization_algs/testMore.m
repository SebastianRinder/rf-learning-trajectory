clear variables;
close all;
set(0, 'DefaultFigureVisible', 'off')

%% seed control
seed = rng;
seed = seed.State(1)
%
seed = 439229998
rng(seed);

%% plots folder
rootPlot = '+static_optimization_algs/plots/';
if(~exist(rootPlot,'dir'))
    mkdir(rootPlot);
end

%% objective function
nbGauss = 25;
muRange = 10;
minSigma = 2;
func = static_optimization_algs.RandomGaussianMixtureFunction(nbGauss, muRange, minSigma);
optimizerInput.fun = func;
%%% plotting 
hnd = figure(1);
func.plot();
funSignature = ['func' num2str(seed) '_' num2str(nbGauss) '_' num2str(muRange) '_' num2str(minSigma)];
fname = [rootPlot 'objective_' funSignature];
hgexport(hnd, [fname '.eps']); %this works better than saveas and print

%% optimization algo
optimizers = {...
%      static_optimization_algs.RepsBandits; ...
    static_optimization_algs.More;
};

all_perfs = {};
all_kls = {};
all_signatures = {};

% init distribution and hyper-params
optimizerInput.initVar = 50;
mu = [0 0]; covC = eye(2) * optimizerInput.initVar;
optimizerInput.initDistrib = static_optimization_algs.Normal(mu, covC);
optimizerInput.epsiKL = .01;
optimizerInput.nbSamplesPerIter = 10;
optimizerInput.nbInitSamples = 10;
optimizerInput.nbIter = 100;
optimizerInput.maxIterReuse = 0;
optimizerInput.useImportanceSampling = true;
optimizerInput.entropyReduction = .1;

% This... (if you want to fix the reularization param)
optimizerInput.regularization = 1e-6;
% ... or this (hyperParam optim)
% optimizerInput.regularization.hyperParamList = 10.^[-6:2];
% optimizerInput.regularization.nbCrossValidation = 50;
% optimizerInput.regularization.validationSetRatio = .15;

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
    [all_perfs{end+1}, all_kls{end+1}] = optimizers{k}.optimizeStruct(optimizerInput);
    toc
end

%% performance plotting
fname = [rootPlot 'perfOn_' funSignature '_' all_signatures{1}];
hnd = figure(5);
hold on;
for k = 1:length(optimizers)
    plot(all_perfs{k}(:, 1), all_perfs{k}(:, 2));
    perfs = all_perfs{k};
end
legHnd = legend(all_signatures{:}, 'Location','southeast');
set(legHnd, 'interpreter', 'none');
set(legHnd, 'FontSize', 7);
hgexport(hnd, [fname '.eps']);

%% kl plotting
fname = [rootPlot 'klOn_' funSignature '_' all_signatures{1}];
hnd = figure(6);
for k = 1:length(optimizers)
    semilogy(all_perfs{k}(:, 1), all_kls{k});
    hold on;
end
hold off;
legHnd = legend(all_signatures{:}, 'Location','southeast');
set(legHnd, 'interpreter', 'none');
set(legHnd, 'FontSize', 7);
hgexport(hnd, [fname '.eps']); 

