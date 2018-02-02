classdef globalBO
    methods (Static)
        function ret = optimze(func, trial)
            allY = [];
            knownY = [];
            hyperTrace = [];
            samples = [];
            trajectories = cell(1,1);

            ub = func.opts.ub;
            lb = func.opts.lb;

            initialSamplesCount = 10;

            bsf = 0;
            for i = 1:initialSamplesCount
                X = randBound(lb, ub, initialSamplesCount);
                if minDist(X,lb,ub) > bsf
                    bsf = minDist(X,lb,ub);
                    samples = X;
                end
            end

            for i=1:initialSamplesCount
                for j = 1:func.opts.trajectoriesPerSample
                    [knownY(i,j), trajectories(i,j)] = func.eval(samples(i,:));
                end
            end

        %     D = func.opts.distanceMat(randBound(lb,ub,10000), samples, trajectories, false, func.opts);
        %     [~,sigmal1] = func.opts.scaleKernel(D,[]);
        %     D = func.opts.distanceMat(samples, samples, trajectories, false, func.opts);
        %     [~,sigmal2] = func.opts.scaleKernel(D,[]);
        % %     randSamples = randBound(lb,ub,300);
        % %     D = func.opts.distanceMat(randSamples, randSamples, trajectories, true, func.opts);
        % %     [~,sigmal] = func.opts.scaleKernel(D,[]);
        %     func.opts.hyper = [0,log(mean([sigmal1,sigmal2]))];

            for i = initialSamplesCount+1:initialSamplesCount+func.opts.bayOptSteps
                y = (knownY - max(knownY))./std(knownY);
        %         y = knownY ./ func.opts.timeSteps;

                D = func.opts.distanceMat(samples, samples, trajectories, false, func.opts);
                if func.opts.hyperOptimize
                    func.opts.hyper = optimizeHyper(samples,y,D, func.opts);
                end
                hyperTrace = [hyperTrace; func.opts.hyper];

                [L, alpha] = getLowerCholesky(D, y, false, func.opts.noiseVariance);

                if func.opts.useMaxMean
                    negGPMean = @(testX) -gaussianProcess(testX, samples, trajectories, L, alpha, func);
                    [~,negMaxMean] = globalMinSearch(negGPMean, lb, ub, func.opts.useGADSToolbox, false);
                    bestY = -negMaxMean;
                else
                    bestY = max(y);
                end           

                negAcqFcn = @(testX) -expectedImprovement(testX, samples, trajectories, L, alpha, func, bestY);
                samples(i,:) = globalMinSearch(negAcqFcn, lb, ub, func.opts.useGADSToolbox, false);
                if func.opts.acquisitionPlot
                    selectFigure('Expected Improvement values (sorted)');
                    xplot = randBound(lb,ub,10000);
                    yplot = negAcqFcn(xplot);
                    plot(sort(-yplot));
                    pause(0.1);
                    hold on;
                end

                for j = 1:func.opts.trajectoriesPerSample
                    [allY(i,j), trajectories(i,j)] = func.eval(samples(i,:));
                end
                knownY(i,1) = mean(allY(i,:));
                if func.opts.verbose == 1
                    disp(['trial ',num2str(trial) ,...
                        ', step ',num2str(i-initialSamplesCount) ,...
                        ', cr ',num2str(knownY(i,1)),...
                        ', sf ',num2str(exp(hyperTrace(end,1))),...
                        ', sl ',num2str(exp(hyperTrace(end,2)))...
                        ]);
                end
            end
            ret = knownY;
        end
    end
end