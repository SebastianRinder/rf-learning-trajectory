clear variables;
close all;

% local bayes vs more 
rootPlot = '+static_optimization_algs/xp_smoothGauss_10speri_20var_min5/';
if(~exist(rootPlot,'dir'))
    mkdir(rootPlot);
end

%%% generating the functions
%some randomly generating seeds
seeds = [3893652211, 642246116, 441759345, 1589975792, 4183876962,...
    685004504, 267927554, 257413324, 3096065520, 3821287115, 1700735318];
% seeds = seeds(1:2)

nbGauss = 25;
muRange = 10;
minSigma = 2;
nbFun = length(seeds);
funs = cell(nbFun, 1);

for k = 1:nbFun
    func = static_optimization_algs.RandomGaussianMixtureFunction(nbGauss, muRange, minSigma, 2, seeds(k));
    funs{k} = func;
    %%% plotting objective
    hnd = figure(1);
    func.plot();
    funSignature = func.getSignature;
    fname = [rootPlot 'objective_' funSignature];
    hgexport(hnd, [fname '.eps']); %this works better than saveas and print
end

%%% the algorithms
algs =  {static_optimization_algs.More; static_optimization_algs.LocalBayes;};

%%% settings
%common settings
nbAlgs = length(algs);
settings = cell(nbAlgs, 1);
initVar = 20;
mu = [0 0];
covC = eye(2) * initVar;

for k = 1:nbAlgs
    settings{k}.initVar = initVar;
    settings{k}.initDistrib = static_optimization_algs.Normal(mu, covC);
    
    settings{k}.epsiKL = .05;
    settings{k}.entropyReduction = .05;
    
    settings{k}.nbSamplesPerIter = 5;
    settings{k}.nbInitSamples = 5;
    settings{k}.nbIter = 100;    
    settings{k}.maxIterReuse = 40;    
end

%more settings
settings{1}.useImportanceSampling = 1;
settings{1}.regularization = 1e-6;

%lb settings
settings{2}.gpStuffPath = 'GPstuff-4.7/';
settings{2}.nbPStarSamples = 100;
settings{2}.minProbaReuse = inf;
settings{2}.gpHyperOption = 'MCMC';
settings{2}.samplingOption = 'Acquisition';
settings{2}.kernelType = 'sexp';
settings{2}.featureFunction = [];
settings{2}.featureName = 'noFeature';
settings{2}.yCenteringType = 'min';

%%%%%%%%%%%%%%%%% calling the optimizers %%%%%%%%%%%%%%%%%%%
save(fullfile(rootPlot, 'xp_settings'), 'funs', 'algs', 'settings');
[all_perfs, all_evals] = static_optimization_algs.batchTestOptim(funs, algs, settings, rootPlot);
save(fullfile(rootPlot, 'allPerfsAndEvals'), 'all_perfs', 'all_evals');

% plottings
alg_signatures = cell(nbAlgs, 1);
for algi = 1:nbAlgs
    alg_signatures{algi} = algs{algi}.getSignature(settings{algi});
end

for k = 1:nbFun
    fname = fullfile(rootPlot, ['perfOn_' funs{k}.getSignature()]);
    hnd = figure;
    hold on;
    for algi = 1:nbAlgs
        plot(all_perfs{k, algi}(:, 1), all_perfs{k, algi}(:, 2));        
    end    
    legHnd = legend(alg_signatures{:}, 'Location','southeast');
    set(legHnd, 'interpreter', 'none');
    set(legHnd, 'FontSize', 7);
    hgexport(hnd, [fname '.eps']);
end
