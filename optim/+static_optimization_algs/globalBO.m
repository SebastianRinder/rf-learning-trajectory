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
            
            global nv;
            nv = func.opts.noiseVariance;
            func.opts.hyper = log([1 1]);
            hypOptIter = 0;
            if strcmp(func.opts.actionSpace,'continuous')
                func.opts.hyperOptimize = 1;
                hyperOptimize = 1;
            else
                func.opts.hyperOptimize = 0;
                hyperOptimize = 0;
            end
            for i = initialSamplesCount+1:initialSamplesCount+func.opts.bayOptSteps
%                 stdy = max(std(knownY),1e-6);
                y = knownY;
%                 y = (knownY - max(knownY))./stdy;
%                 y = (knownY - min(knownY))./stdy;
%                 y = (knownY - mean(knownY))./stdy;
        %         y = knownY ./ func.opts.timeSteps;

                D = func.opts.distanceMat(samples, samples, trajectories, 'kk', func.opts); %known known
                if hyperOptimize && func.opts.hyperOptimize
                    func.opts.hyper = optimizeHyper(samples, y, D, func.opts);
%                     func.opts.hyper(2) = log(1400);
                    hyperOptimize = 0;
                end
                hyperTrace = [hyperTrace; func.opts.hyper];
                
                K = func.opts.scaleKernel(D, func.opts.hyper);
                [L, alpha] = getLowerCholesky(K, y, false, func.opts.noiseVariance);

                if func.opts.useMaxMean
                    negGPMean = @(testX) -gaussianProcess(testX, samples, trajectories, L, alpha, func);
                    [~,negMaxMean] = globalMinSearch(negGPMean, lb, ub, func.opts.useGADSToolbox, false);
                    bestY = -negMaxMean;
                else
                    bestY = max(y);
                end           

                negAcqFcn = @(testX) -expectedImprovement(testX, samples, trajectories, L, alpha, func, bestY);
%                 tic
                [samples(i,:), negEI] = globalMinSearch(negAcqFcn, lb, ub, func.opts.useGADSToolbox, false);
                if -negEI <= 1e-10
                    warning('off','backtrace');
                    warning('optimum of EI near 0');
                    warning('on','backtrace');
                end
                hypOptIter = hypOptIter + 1;
                if -negEI <= 1e-6 || hypOptIter == 10
                    hypOptIter = 0;
                    hyperOptimize = 1;
                end

                if func.opts.acquisitionPlot
                    selectFigure('random Expected Improvement values (sorted)');
                    xplot = randBound(lb,ub,10000);
                    yplot = negAcqFcn(xplot);
                    plot(sort(-yplot));
                    if min(yplot) < negEI
                        warning('off','backtrace');
                        warning('optimum of EI worse than best of random 10000');
                        warning('on','backtrace');
                    end
                    pause(0.1);
%                     hold on;
%                     if sum(yplot==0) > 9000 || -negEI <= 1e-10
%                         func.opts.hyperOptimize = 1;
%                     end
                end

                for j = 1:func.opts.trajectoriesPerSample
                    [allY(i,j), trajectories(i,j)] = func.eval(samples(i,:));
                end
                knownY(i,1) = mean(allY(i,:));
                if func.opts.verbose == 1
                    disp(['trial ',num2str(trial) ,...
                        ', step ',num2str(i-initialSamplesCount) ,...
                        ', cr ',num2str(knownY(i,1)),...
                        ', sf ',num2str(exp(hyperTrace(end,1)))...
                        ', sl ',num2str(exp(hyperTrace(end,2)))...
                        ]);
                end
            end
            ret = knownY;
        end
    end
end