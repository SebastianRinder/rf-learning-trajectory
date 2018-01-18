function [experiment] = changeExperimentPath(experiment2, path)
    experiment = experiment2;
    if(path(end) ~= '/')
        path(end+1) = '/';
    end
    experiment.path = path;
    [~, setting] = fileparts(experiment.experimentPath);
    experiment.experimentPath = [path, setting];
    if(setting(end) ~= '/')
        setting(end+1) = '/';
    end
    for i = 1:length(experiment.evaluations)
        [~, eval] = fileparts(experiment.evaluations{i}.path);
        experiment.evaluations{i}.path =  [path, setting, eval];
    end
end
