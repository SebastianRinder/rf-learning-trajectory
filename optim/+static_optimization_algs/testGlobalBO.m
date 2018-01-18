%% author: Sebastian Rinder

function ret = testGlobalBO()
    isCluster = 1; %0 if debugging
    if ~isCluster
        ret = licenceLoopGlobalBO(0);
        return;
    else
        loop = true;
        cores = 32;
        while loop
            try
                delete(gcp('nocreate')); %shut down previously created parpool
                parpool(cores);
                ret = licenceLoopGlobalBO(cores);
                loop = false;
            catch me
                if ~strcmp(me.identifier,'parallel:cluster:LicenseUnavailable')
                    hrs = datestr(now,'_dd-mm-yyyy_HH-MM');
                    save(['results/error',hrs,'.mat'],'me');
                    disp(me.message);
                    loop = false;
                else
                    pause(30);
                end
            end
        end
    end
end

function ret = licenceLoopGlobalBO(cores)
    %dbstop if error

    kernel = {'sexp','matern52','trajectory'};
    
%     platform = 'pygym';
    platform = 'matlab';
    
    env = 'cartPole';
    
    addpath('auxiliary');
    addpath('cartPole');
%     addpath('openAIGym');

    trials = 32;
    bayOptSteps = 200;
    initialSamplesCount = 10;

%     ub = 1.*ones(1,func.opts.dim);
%     lb = -1.*ub;
%     func = problemOptions(kernel{kernelIdx},platform,env);    
%     noiseVariance = getNoiseVariance(func,lb,ub);
    
%     noiseVariance = 0.2;
    noiseVariance = 1e-6;
%     noiseVariance = 0.1;
    nnn = logspace(-8,3.5,5);

    useMaxMean = 0; % if the environment is very noisy use maximum mean of the GP instead of max(knownY) in the expected improvement (Brochu 2010)
    
    ret = cell(trials,1);
%     for kernelIdx = 1:3
    kernelIdx = 3;
    for noiseIdx = 1:length(nnn)
        noiseVariance = nnn(noiseIdx);
        parfor (trial = 1:trials, cores)
%         for trial = 1:trials
            tic
            allY = [];
            knownY = [];
            hyperTrace = [];
            samples = [];
            trajectories = cell(1,1);

            func = problemOptions(kernel{kernelIdx},platform,env);
            func.opts.hyper = [0,0];
            func.opts.trajectoriesPerSample = 1;
            func.opts.hyperOptimize = 0;
            func.opts.hyperPlot = 0;
            func.opts.acquisitionPlot = 0;
            func.opts.useGADSToolbox = 0;
            func.opts.sigmaNoiseSquared = noiseVariance;
            ub = 1.*ones(1,func.opts.dim);
            lb = -1.*ub;

            bsf = 0;
            for i = 1:100
                X = randBound(lb, ub, initialSamplesCount);
                if minDist(X,lb,ub) > bsf
                    bsf = minDist(X,lb,ub);
                    samples = X;
                end
            end

            for i=1:initialSamplesCount
                for j = 1:func.opts.trajectoriesPerSample
                    [allY(i,j), trajectories(i,j)] = func.eval(samples(i,:));
                end
                knownY(i,1) = mean(allY(i,:));
            end

            D = func.opts.distanceMat(randBound(lb,ub,10000), samples, trajectories, false, func.opts);
            [~,sigmal] = func.opts.scaleKernel(D,[]);
            func.opts.hyperLb = [-3,-20];
            func.opts.hyperUb = [3,20];
    %         func.opts.hyperLb = log(sigmal)-10;
    %         func.opts.hyperUb = log(sigmal)+10;
            func.opts.hyper = [0,log(sigmal)];

            for i = initialSamplesCount+1:initialSamplesCount+bayOptSteps
                y = (knownY - mean(knownY))./std(knownY);
                D = func.opts.distanceMat(samples, samples, trajectories, false, func.opts);
                if func.opts.hyperOptimize
                    func.opts.hyper = optimizeHyper(samples,y,D, func.opts);
                end
                hyperTrace = [hyperTrace; func.opts.hyper];

                [L, alpha] = getLowerCholesky(D, y, false, func.opts.sigmaNoiseSquared);

                if useMaxMean
                    negGPMean = @(testX) -gaussianProcess(testX, samples, trajectories, L, alpha, func);
                    [~,negMaxMean] = globalMinSearch(negGPMean, lb, ub, func.opts.useGADSToolbox, false);
                    bestY = -negMaxMean;
                else
                    bestY = max(y);
                end           

    %             func.opts.sigmaNoiseSquared = 0.2;%1e-6;
                negAcqFcn = @(testX) -expectedImprovement(testX, samples, trajectories, L, alpha, func, bestY);
                samples(i,:) = globalMinSearch(negAcqFcn, lb, ub, func.opts.useGADSToolbox, false);
                if func.opts.acquisitionPlot
                    selectFigure('Expected Improvement values (sorted)');
                    xplot = randBound(lb,ub,10000);
                    yplot = negAcqFcn(xplot);
                    plot(sort(-yplot));
                    pause(0.1);
                end

                for j = 1:func.opts.trajectoriesPerSample
                    [allY(i,j), trajectories(i,j)] = func.eval(samples(i,:));
                end

                knownY(i,1) = mean(allY(i,:));
                disp(['krnl ',num2str(kernel{kernelIdx}) ,...
                    ', trial ',num2str(trial) ,...
                    ', step ',num2str(i-initialSamplesCount) ,...
                    ', cr ',num2str(knownY(i,1)),...
                    ', sf ',num2str(exp(hyperTrace(end,1))),...
                    ', sl ',num2str(exp(hyperTrace(end,2)))...
                    ]);

            end

            %             ret{trial,1}.samples = samples;
            %             ret{trial,1}.trajectories = trajectories;
            ret{trial,1}.hyperTrace = hyperTrace;
            ret{trial,1}.knownY = knownY;
            ret{trial,1}.timeTakenSeconds = toc;
            disp(num2str(toc/60));
        end
        hrs = datestr(now,'dd-mm-yyyy_HH-MM');
        saveStr = sprintf('results/%s_%s_%s_%0.0e_%s.mat', env, platform, kernel{kernelIdx}, noiseVariance, hrs);
        save(saveStr,'ret');
    end
    
    delete(gcp('nocreate'));
end
