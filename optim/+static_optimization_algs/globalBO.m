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

            for i = initialSamplesCount+1:initialSamplesCount+func.opts.bayOptSteps
                stdy = max(std(knownY),1e-6);
                y = (knownY - max(knownY))./stdy;
                y = (knownY - mean(knownY))./stdy;
        %         y = knownY ./ func.opts.timeSteps;

                D = func.opts.distanceMat(samples, samples, trajectories, 'kk', func.opts); %known known
%                 if func.opts.hyperOptimize
%                     func.opts.hyper = optimizeHyper(samples, y, D, func.opts);
%                 end
                hyperTrace = [hyperTrace; func.opts.hyperl.uk, func.opts.hyperl.kk, func.opts.hyperl.uu];
                
                K = func.opts.scaleKernel(D, [0, func.opts.hyperl.kk]);
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
                samples(i,:) = globalMinSearch(negAcqFcn, lb, ub, func.opts.useGADSToolbox, false);
%                 [xxx,yyy] = globalMinSearch(negAcqFcn, lb, ub, func.opts.useGADSToolbox, false)
%                 toc
%                 tic
%                 [xxx,yyy] = globalMinSearch(negAcqFcn, lb, ub, ~func.opts.useGADSToolbox, false)
%                 toc
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
                        ', sluk ',num2str(exp(hyperTrace(end,1)))...
                        ', slkk ',num2str(exp(hyperTrace(end,2)))...
                        ', sluu ',num2str(exp(hyperTrace(end,3)))...
                        ]);
                end
            end
            ret = knownY;
        end
    end
end