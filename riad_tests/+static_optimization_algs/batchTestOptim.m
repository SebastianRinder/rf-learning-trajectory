%run all alg{i} with setting{i} on all functions funs{j}
function [all_perfs, all_evals] = batchTestOptim(funs, algs, settings, rootPlot)

all_perfs = cell(length(funs), length(algs));
all_evals = cell(length(funs), length(algs));

for fi = 1:length(funs)
    seed = rng;
    seedStartOpt = seed.State(2); % why two?
    for algi = 1:length(algs)
        rng(seedStartOpt);
        if(exist('rootPlot', 'var'))
            settings{algi}.videoFile = VideoWriter(fullfile(rootPlot, [algs{algi}.getSignature(settings{algi}) '_' funs{fi}.getSignature]));
            settings{algi}.videoFile.FrameRate = 4;
        end
        [all_perfs{fi, algi}, all_evals{fi, algi}] = algs{algi}.optimizeStruct(settings{algi}, funs{fi});
    end
end

end

