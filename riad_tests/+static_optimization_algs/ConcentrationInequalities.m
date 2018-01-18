classdef ConcentrationInequalities
    methods (Static)
        function [lb, hatWR] = andersonLowerBound(vals, logImportanceProbas, logProbas, conf)
            z = sort(vals .* exp(logProbas - logImportanceProbas));
            hatWR = mean(z);
            cst = sqrt(log(2/conf) / (2 * length(z)));
            weightDiff = min(((1:length(z))' - 1) ./ length(z) + cst, 1);
            lb = z(end) - (z - [0; z(1:end-1)])' * weightDiff;
        end
        
        function [lb, hatWR, sigWR] = studentLowerBound(vals, logImportanceProbas, logProbas, conf)
            weightedR = vals .* exp(logProbas - logImportanceProbas);
            hatWR = mean(weightedR);
            sigWR = std(weightedR);
            lb = hatWR - sigWR * tinv(1-conf, length(vals) - 1) / sqrt(length(vals));
        end
        
        function [selection, tooLess] = trimImprobableSampels(logProbas, logImportanceProbas, ratio, minSamples)
            [sWeights, is] = sort(exp(logProbas - logImportanceProbas), 'descend');
            firstId = find(cumsum(sWeights / sum(sWeights)) > ratio, 1, 'first');
            if(firstId < minSamples + 1)
                tooLess = 1;
                selection = [];
            else
                tooLess = 0;
                selection = is(1:firstId-1);
            end
        end
        
        function [weightedRewards, tooLess] = trimImprobable(logProbas, logImportanceProbas, vals, ratio, minSamples)
            [sWeights, is] = sort(exp(logProbas - logImportanceProbas), 'descend');
            firstId = find(cumsum(sWeights / sum(sWeights)) > ratio, 1, 'first');
            tooLess = 0;
            if(firstId < minSamples + 1)
                weightedRewards = [];
                tooLess = 1;
                return;
            else
                is = is(1:firstId-1);
                logImportanceProbas = logImportanceProbas(is);
                maxLogProbas = max(logImportanceProbas);
                logImportanceProbas = logImportanceProbas - (log(sum(exp(logImportanceProbas - maxLogProbas))) + maxLogProbas);
                sWeights = exp(logProbas(is) - logImportanceProbas);
                weightedRewards = vals(is) .* sWeights;
            end
        end        
    end
end
