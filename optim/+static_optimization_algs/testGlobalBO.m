%% author: Sebastian Rinder


function ret = testGlobalBO()
    loop = 1;
    while loop
        try
            delete(gcp('nocreate')); %shut down previously created parpool
            parpool(32);
            ret = licenceLoopGlobalBO();
            loop = 0;
        catch me
            %disp(me.message);
            pause(30);
            loop = 1;
        end
    end
end

function ret = licenceLoopGlobalBO()
    %dbstop if error
    %xp = 'CartpoleSexp';
    xp = 'CartpoleTraj';
    
    delete(['ret',xp,'.mat']);

    addpath('auxiliary');
    addpath('cartPole');
               
    trials = 32;
    bayOptSteps = 200;
    initialSamples = 10;
        
    useMaxMean = 0; % if the environment is very noisy use maximum mean of the GP instead of max(knownY) in the expected improvement (Brochu 2010)
    
    ret = cell(trials,1);
    
    parfor trial = 1:trials
        tic
        allY = [];
        knownY = [];
        hyperTrace = [];
        samples = [];
        trajectories = cell(1,1);
        
        func = problemOptions('trajectory','cartPole');
        func.opts.hyper = [0,0];
        func.opts.trajectoriesPerSample = 1;
        func.opts.hyperPlot = 0;
        func.opts.acquisitionPlot = 0;
        func.opts.useGADSToolbox = 0;
        
        ub = 1.*ones(1,func.opts.dim);
    	lb = -1.*ub;
        
        for i=1:initialSamples
            samples(i,:) = randBound(lb, ub, 1);
            for j = 1:func.opts.trajectoriesPerSample
                [allY(i,j), trajectories(i,j)] = func.eval(samples(i,:));
            end
            knownY(i,1) = mean(allY(i,:));
        end

        for i = initialSamples+1:initialSamples+bayOptSteps
            y = (knownY - mean(knownY))./std(knownY);
            %func.opts.sigmaNoiseSquared = max(mean(std(allY,0,2)).^2, 1e-6);
            %func.opts.sigmaNoiseSquared = max(std(knownY)/sqrt(2), 1); %from matlab
            func.opts.sigmaNoiseSquared = 0.2;  %from matlab bayesopt: std(y) / 5
            D = func.opts.distanceMat(samples, samples, trajectories, false, func.opts);
            func.opts.hyper = optimizeHyper(samples,y,D, func.opts);
%             func.opts.hyper = [0,0];
            hyperTrace = [hyperTrace; func.opts.hyper];
            
            [L, alpha] = getLowerCholesky(D, y, false, func.opts.sigmaNoiseSquared);
            
            if useMaxMean 
                negGPMean = @(testX) -gaussianProcess(testX, samples, trajectories, L, alpha, func);
                [~,negMaxMean] = globalMinSearch(negGPMean, lb, ub, func.opts.useGADSToolbox, false);
                bestY = -negMaxMean;
            else
                bestY = max(y);
            end           
            
            func.opts.sigmaNoiseSquared = 0.2;%1e-6;
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
            disp(['trial : ',num2str(trial) ,...
                ' | step : ',num2str(i-initialSamples) ,...
                ' | cr: ',num2str(knownY(i,1)),...
                ' | sigmaf: ',num2str(exp(hyperTrace(end,1))),...
                ' | sigmal: ',num2str(exp(hyperTrace(end,2)))...
                ]);

%             if ~isequal(func.opts.visualize, 'none')
%                 [~,bestTrajectoryIdx] = max(allY(i,:));
%                 [~,bestObjectiveIdx] = max(knownY);
%                 if isequal(func.opts.visualize, 'best') && i == bestObjectiveIdx ||...
%                         isequal(func.opts.visualize, 'all')
%                     traj = trajectories(i,bestTrajectoryIdx);
%                     opts.visFcn(traj,opts.bounds);
%                 end
%             end
%             fprintf(' %d',i-initialSamples);
        end
        
        %             ret{trial,1}.samples = samples;
        %             ret{trial,1}.trajectories = trajectories;
        %             ret{trial,1}.hyperTrace = hyperTrace;
        ret{trial,1}.knownY = knownY;
        ret{trial,1}.timeTakenSeconds = toc;
        disp(num2str(toc/60));
    end
    save(['ret',xp,'.mat'],'ret');
    
    delete(gcp('nocreate'));
end
