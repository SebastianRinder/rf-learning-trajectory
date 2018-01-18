clear all;
pre = '/home/akrour/postdoc/code/policysearchtoolbox/+Experiments/data/test/TableTennis_TrajectoryBased_SimpleAlg/Table_Tennis_201507242016_01/eval001/trial001/iter_';
post = '.mat';
reward = Environments.SL.Tasks.SLTableTennisReward;
maxIter = 400000;
topRewards = zeros(maxIter, 1);
argMaxes = zeros(maxIter, 1);
idx = 1;
for n = 1:10:maxIter
    lower = int2str(n);
    while(length(lower) < 5)
        lower = ['0' lower];
    end
    higher = int2str(n+9);
    while(length(higher) < 5)
        higher = ['0' higher];
    end
    fname = [pre lower '_' higher post];
    disp(fname);
    clear iter*;
    allvars = load(fname);
    for f = fieldnames(allvars)'        
        if(isfield(allvars.(f{1}), 'data'))
            topRewards(idx) = reward.getTrajectoryReward(allvars.(f{1}).data.steps.SLstates);
            currentIter = str2num(f{1}(1, 5:end));
            argMaxes(idx) = currentIter;
            idx = idx + 1;   
        end
    end    
end
topRewards = topRewards(1:idx-1);
argMaxes = argMaxes(1:idx-1);
[topRewards, indxs] = sort(topRewards, 'descend');
argMaxes = argMaxes(indxs);
save(['readIterResults' datestr(now)], 'topRewards', 'argMaxes');
