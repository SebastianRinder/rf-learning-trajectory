classdef RepsTrueGrid
    methods (Static)
        function sign = getSignature(optimizerInput)
            sign = ['NonParamRepsGrid' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbBinPerDim) '_' num2str(optimizerInput.initVar)];
        end
        
        function [perf] = optimizeStruct(optimizerInput)
            perf = static_optimization_algs.RepsTrueGrid.optimize(optimizerInput.initDistrib, optimizerInput.epsiKL,...
                optimizerInput.nbIter, optimizerInput.nbBinPerDim, optimizerInput.nbSamplesPerIter, ...
                optimizerInput.fun, optimizerInput.videoFile);
        end
       
        function [perf] = optimize(initDistrib, epsiKL, nbIter, nbBinPerDim, nbSamplesPerIter, fun, videoFile)
            perf = zeros(nbIter, 2);
            % won't bother with dimensions higher than 2
            [limitInf, limitSup] = fun.getRange();
            line1 = limitInf(1):(limitSup(1)-limitInf(1))/(nbBinPerDim-1):limitSup(1);
            line2 = limitInf(2):(limitSup(2)-limitInf(2))/(nbBinPerDim-1):limitSup(2);
            [x1, x2] = meshgrid(line1, line2);
            samples = [x1(:) x2(:)];
            sampleInitProbas = exp(initDistrib.getLogProbas(samples));
            vals = fun.eval(samples);
            diffVals = vals - max(vals);
            sum1OverEta = 0;
            %% animation
            open(videoFile);
            
            for iter = 1:nbIter
                %% evaluation
                sampleWeights = exp(diffVals * sum1OverEta);
                sampleWeights(isnan(sampleWeights)) = 1;
                sampleWeights = sampleWeights .* sampleInitProbas;
                sampleWeights = sampleWeights / sum(sampleWeights);
                
                newPerf = vals' * sampleWeights / sum(sampleWeights);
                perf(iter, :) = [((iter-1) * nbSamplesPerIter) newPerf];               
               
                %% current distrib plotting
                frame = static_optimization_algs.RepsTrueGrid.getFrame(sampleWeights, line1, line2);
                writeVideo(videoFile, frame);

                
                %% update
                eta = static_optimization_algs.RepsBandits.optimizeDual(vals, epsiKL, sampleWeights);
                sum1OverEta = sum1OverEta + 1/eta;
            end            
            close(videoFile);
        end
        
        function frame = getFrame(sampleWeights, line1, line2)
            currPlot = figure;
            contour(line1, line2, reshape(sampleWeights, length(line2), length(line1)));
            frame = getframe(gcf,[0 0 560 420]);
            close(currPlot);
        end
        
    end
end

