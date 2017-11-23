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
clear variables;
%clc;
% set(0, 'DefaultFigureVisible', 'off')

%% seed control
rng('default');
rng('shuffle');
seed = rng;
seed = seed.State(1)
% seed = 526950115
% seed = 142671073;
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
all_perfs = {};
all_signatures = {};

% init distribution and co.
optimizerInput.fun = func;
optimizerInput.initVar = 10;
mu = [0 0]; covC = eye(2) * optimizerInput.initVar;
optimizerInput.xinit = 'varargin{1}{2}.opts.initDistrib.getSamples(1)';
% optimizerInput.xinit = [0; 0];
optimizerInput.initDistrib = static_optimization_algs.Normal(mu, covC);
optimizerInput.maxEvals = 2000; %cmaes wrapper only depends on this.
optimizerInput.nbSamplesPerIter = 6;
optimizerInput.maxIterReuse = 30;
optimizerInput.maxIter = 20;

% opt_alg = static_optimization_algs.CMAESWrapper(optimizerInput, func);
% opt_alg = static_optimization_algs.CMAESGPMean(optimizerInput, func);
opt_alg = static_optimization_algs.CMAESThompson(optimizerInput, func);

seed = rng;
seedStartOpt = seed.State(2)
rng(seedStartOpt);
all_signatures{end+1} = [opt_alg.getSignature() '_' funSignature];

videoFile = VideoWriter([rootPlot 'policy_search_' all_signatures{end}]);%, 'Uncompressed AVI');
videoFile.FrameRate = 4;
opt_alg.video = videoFile;

disp(['starting ' all_signatures{end}]);
tic
opt_alg.optimize()
[all_perfs{end+1}] = opt_alg.allEvals;
toc

%% performance plotting
fname = [rootPlot 'perfOn_' funSignature '_' all_signatures{1}];
hnd = figure(5);
hold on;
plot(all_perfs{1}(:, 1));
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

