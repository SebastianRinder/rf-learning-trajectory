classdef BayesOptLib
    % wraper to the bayesopt library
    
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['BayesOptLib' '_'  num2str(optimizerInput.nbInitSamplesBO) '_' ...
                num2str(optimizerInput.nbEvalsBO) '_' num2str(optimizerInput.maxNbSamplesBO)];
        end
        
        %%
        function [allVals] = optimizeStruct(optimizerInput, func)
            if(~isfield(optimizerInput, 'videoFile'))
                optimizerInput.videoFile = [];
            end
            
            [allVals] = static_optimization_algs.BayesOptLib.optimize(optimizerInput.functionBounds, optimizerInput.nbInitSamplesBO,...
                optimizerInput.nbEvalsBO, optimizerInput.maxNbSamplesBO, func, optimizerInput.videoFile);
        end
        
        %% main function
        function [allVals] = optimize(functionBounds, nbInitEvals, nbTotalEvals, maxNbSamples, fun, videoFile)
            %%% bayes_optim
            params.n_iterations = nbTotalEvals - nbInitEvals;
            params.n_init_samples = nbInitEvals;
            params.crit_name = 'cEI';
            params.surr_name = 'sGaussianProcessML';
            params.noise = 1e-6;
            params.kernel_name = 'kSEARD';
            params.verbose_level = 1;
            params.log_filename = 'matbopt.log';
            params.l_type = 'L_MCMC';
            params.l_all = 1;

            disp('Starting bayes optim. Might take a while');
            
            tic
            bayesoptcont(@(x)-fun.eval(x), size(functionBounds, 2), params, functionBounds(1, :), functionBounds(2, :));
            toc

            allSamples = fun.allSamples;
            allVals = fun.allEvals;
            
            %%% animation
            if(~isempty(videoFile))
                open(videoFile);
                for iter = 1:nbTotalEvals
                    %%% current distrib plotting
                    if(~isempty(videoFile))
                        frame = static_optimization_algs.BayesOptLib.getFrame(fun, allSamples(1:iter, :), iter);
                        writeVideo(videoFile, frame);
                    end
                end                
                close(videoFile);
            end            
        end
        
        function frame = getFrame(fun, samples, iter)
            persistent currPlot;
            if(isempty(currPlot))
                currPlot = figure;
            end
            
            %function and distrib
            fun.plot();
            hold on;
            if(~isempty(samples))
                plot(samples(:, 1), samples(:, 2), '*r');
            end
            hold off;  
            text(-5, -5, ['iteration = ' num2str(iter)]);
            frame = getframe(gcf,[0 0 560 420]);

        end
    end
end
