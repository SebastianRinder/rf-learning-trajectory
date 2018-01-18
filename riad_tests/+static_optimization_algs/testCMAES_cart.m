close all;
clear variables;
%% seed control
rng('default');
rng('shuffle');
seed = rng;
seed = seed.State(1)
seed = 526950115
rng(seed);
rootPlot = '+static_optimization_algs/plots/';
if(~exist(rootPlot,'dir'))
    mkdir(rootPlot);
end

%% objective function
addpath('cartPole')
func = cartPole();
funSignature = ['cartpole'];
fname = [rootPlot 'objective_' funSignature];

%% optimization algo
all_perfs = {};
all_signatures = {};

% init distribution and co.
optimizerInput.fun = func;
optimizerInput.initVar = 1;
mu = zeros(1, func.opts.dim);
covC = eye(func.opts.dim) * optimizerInput.initVar;
% optimizerInput.xinit = 'varargin{1}{2}.opts.initDistrib.getSamples(1)';
optimizerInput.xinit = mu;
optimizerInput.initDistrib = static_optimization_algs.Normal(mu, covC);
optimizerInput.maxEvals = 2000; %cmaes wrapper only depends on this.
optimizerInput.maxIterReuse = 30;
optimizerInput.maxIter = 100;

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
opt_alg.video = [];

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

