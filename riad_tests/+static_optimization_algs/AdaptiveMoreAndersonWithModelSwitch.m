classdef AdaptiveMoreAndersonWithModelSwitch
    methods (Static)
        %%
        function sign = getSignature(optimizerInput)
            sign = ['AdaptiveMoreAndersonWithModelSwitch' '_' num2str(optimizerInput.epsiKL) '_' num2str(optimizerInput.nbIter) '_' ...
                num2str(optimizerInput.nbSamplesPerIter) '_' num2str(optimizerInput.nbInitSamples) '_' ...
                num2str(optimizerInput.maxIterReuse) '_' num2str(optimizerInput.initVar)];
        end
        %%
        function [perf, kls] = optimizeStruct(optimizerInput)
            avrRewardLB = @(vals, logImportanceProbas, logProbas) static_optimization_algs.ConcentrationInequalities.andersonLowerBound(vals, logImportanceProbas, logProbas, optimizerInput.confidenceStudent);
            [perf, kls] = static_optimization_algs.AdaptiveMoreStudent.optimize(optimizerInput.initDistrib, optimizerInput.minEpsiKL,...
                optimizerInput.nbSamplesPerIter, optimizerInput.nbInitSamples, optimizerInput.nbIter, optimizerInput.maxIterReuse, avrRewardLB, optimizerInput.regularization, optimizerInput.fun, optimizerInput.videoFile);
        end
    end        
end
