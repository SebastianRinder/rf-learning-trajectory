clear variables;
close all;

% local bayes vs more 
rootPlot = '+static_optimization_algs/lb_quad/';
if(~exist(rootPlot,'dir'))
    mkdir(rootPlot);
end

%%% generating the functions
%some randomly generating seeds
seeds = [3893652211, 642246116, 441759345, 1589975792, 4183876962,...
    685004504, 267927554, 257413324, 3096065520, 3821287115, 1700735318];

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
algs =  {static_optimization_algs.LocalBayes;};

%%% settings
%common settings
nbAlgs = length(algs);
settings = cell(nbAlgs, 1);
initVar = 10;
mu = [0 0];
covC = eye(2) * initVar;

for k = 1:nbAlgs
    settings{k}.initVar = initVar;
    settings{k}.initDistrib = static_optimization_algs.Normal(mu, covC);
    
    settings{k}.epsiKL = .05;
    settings{k}.entropyReduction = .05;
    
    settings{k}.nbSamplesPerIter = 10;
    settings{k}.nbInitSamples = 10;
    settings{k}.nbIter = 100;    
    settings{k}.maxIterReuse = 30;    
end

%more settings
% settings{1}.useImportanceSampling = 1;
% settings{1}.regularization = 1e-6;

%lb settings
settings{1}.gpStuffPath = 'GPstuff-4.7/';
settings{1}.nbPStarSamples = 100;
settings{1}.minProbaReuse = 1 - 1e-4;
settings{1}.gpHyperOption = 'MCMC';
settings{1}.samplingOption = 'Acquisition';
settings{1}.kernelType = 'linear';
%settings{1}.featureFunction = @(x) [ones(size(x, 1), 1) x];
settings{1}.featureFunction = @(data) static_optimization_algs.Quadratic.getQuadFeatures(data, 1, 1);
settings{1}.featureName = 'quadratic';

%%%%%%%%%%%%%%%%% calling the optimizers %%%%%%%%%%%%%%%%%%%
save(fullfile(rootPlot, 'xp_settings'), 'funs', 'algs', 'settings');
all_perfs = static_optimization_algs.batchTestOptim(funs, algs, settings, rootPlot);
save(fullfile(rootPlot, 'all_perfs'), 'all_perfs');

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
