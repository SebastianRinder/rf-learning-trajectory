%% author: Sebastian Rinder

function ret = testGlobalBO(cores)
%     isCluster = logical(cores); %0 if debugging
%     if ~isCluster
%         ret = licenceLoopGlobalBO(0);
%         return;
%     else
%         loop = true;
% %         cores = 32;
%         while loop
%             try
%                 delete(gcp('nocreate')); %shut down previously created parpool
%                 parpool(cores);
%                 ret = licenceLoopGlobalBO(cores);
%                 loop = false;
%             catch me
%                 if ~strcmp(me.identifier,'parallel:cluster:LicenseUnavailable')
%                     hrs = datestr(now,'_dd-mm-yyyy_HH-MM');
%                     save(['results/error',hrs,'.mat'],'me');
%                     disp(me.message);
%                     loop = false;
%                 else
%                     pause(30);
%                 end
%             end
%         end
%     end
% end

isCluster =1; %0 if debugging
kernel = {'sexp','matern52','trajectory'};

% platform = 'pygym';
platform = 'matlab';

env = 'cartPole';

addpath('pygym');
addpath('auxiliary');
addpath('cartPole');

for kernelIdx = 3:3
    func = problemOptions(kernel{kernelIdx},platform,env);
    func.opts.kernel = kernel{kernelIdx};
    func.opts.hyper = [0,0];
    func.opts.trajectoriesPerSample = 1;
    func.opts.hyperOptimize = 0;
    func.opts.hyperPlot = 0;
    func.opts.acquisitionPlot = 0;
    func.opts.useGADSToolbox = 0;
    func.opts.noiseVariance = 1e-6;
    func.opts.ub = 1.*ones(1,func.opts.dim);
    func.opts.lb = -1.*func.opts.ub;
    func.opts.bayOptSteps = 200;
    func.opts.initialSamplesCount = 10;
    func.opts.useMaxMean = 0; % if the environment is very noisy use maximum mean of the GP instead of max(knownY) in the expected improvement (Brochu 2010)

    if isCluster
        trials = 32;
        loop = true;
        
        while loop
            try
                delete(gcp('nocreate')); %shut down previously created parpool
                parpool(cores);
                parfor (trial=1:trials, cores)
                    tic;
                    ret{trial,1}.knownY = static_optimization_algs.globalBO(func, trial);
                    ret{trial,1}.timeTakenSeconds = toc;
                    ret{trial,1}.hyper = func.opts.hyper;
                    ret{trial,1}.noiseVariance = func.opts.noiseVariance;
                end
                delete(gcp('nocreate'));        
                loop = false;
            catch me
                if ~strcmp(me.identifier,'parallel:cluster:LicenseUnavailable')
                    hrs = datestr(now,'_dd-mm-yyyy_HH-MM');
                    save(['results/error',hrs,'.mat'],'me');
                    disp(getReport(me))
                    loop = false;
                else
                    pause(3);
%                     for trial=1:trials
%                         tic;
%                         ret{trial,1}.knownY = static_optimization_algs.globalBO(func, trial);
%                         ret{trial,1}.timeTakenSeconds = toc;
%                         ret{trial,1}.hyper = func.opts.hyper;
%                         ret{trial,1}.noiseVariance = func.opts.noiseVariance;
%                     end
%                     loop = false;
                end
            end
        end
    else
        trials = 2;
        for trial=1:trials
            tic;
            ret{trial,1}.knownY = static_optimization_algs.globalBO(func, trial);
            ret{trial,1}.timeTakenSeconds = toc;
            ret{trial,1}.hyper = func.opts.hyper;
            ret{trial,1}.noiseVariance = func.opts.noiseVariance;
        end
    end
    hrs = datestr(now,'dd-mm-yyyy_HH-MM');
    saveStr = sprintf('results/%s_%s_%s_%s_%0.0e_%s.mat', 'global_ei', env, platform, kernel{kernelIdx}, func.opts.noiseVariance, hrs);
    save(saveStr,'ret');
end