% testing the mixture model
clear all;
close all;
clc;
dbstop if error;

% dataManager, has only one entry
dataManager = Data.DataManager('data');
dataManager.addDataEntry('x', 2);
dataManager.finalizeDataManager();

% GaussianLinear in features. Takes no input
simpleGauss = Distributions.Gaussian.GaussianLinearInFeatures(dataManager, 'x', {}, 'bias');
simpleGauss.initObject(); % important (why not call in constructor?)

% samples from first gaussian
simpleGauss.setBias([0 0]);
simpleGauss.setSigma(0.4 * eye(2));
sample{2} = simpleGauss.sampleFromDistribution(10000); 

% samples from second gaussian
simpleGauss.setBias([4 4]);
simpleGauss.setSigma(eye(2));
sample{3} = simpleGauss.sampleFromDistribution(10000); 


hold on;
title ('Plotting samples from two gaussians');
nbBins = 11;
sample{1} = [];

for i=2:length(sample)
    sample{1} = [sample{1}; sample{i}];
end
[n,c] = hist3(sample{1}, [nbBins nbBins]);
contour(c{1}, c{2}, n);


% now sampling from a mixture model
figure;
hold on;

dataManager.addDataEntry('o', 1, 0, 1);
dataManager.finalizeDataManager();
gating = Distributions.Discrete.ConstantDiscreteDistribution(dataManager, 'o', 'gating');
gating.setItemProb([.5 .5]);
settings = Common.Settings();
settings.setProperty('numOptions', 2);
mixtureModel = Distributions.MixtureModel.MixtureModel(dataManager, gating, ...
    @Distributions.Gaussian.GaussianLinearInFeatures, 'x', {}, 'opt');
mixtureModel.initObject();
mixtureModel.getOption(1).setBias([0 0]);
mixtureModel.getOption(1).setSigma(0.4 * eye(2));
mixtureModel.getOption(2).setBias([4 4]);
mixtureModel.getOption(2).setSigma(eye(2));

sampleMix{2} = mixtureModel.sampleFromDistribution(20000);

sampleMix{1} = [];
title ('Plotting samples from mixture of two gaussians');
for i=2:length(sampleMix)
    sampleMix{1} = [sampleMix{1}; sampleMix{i}];
end

% hist3(sampleMix{i});
[n,c] = hist3(sampleMix{i}, [nbBins nbBins]);
contour(c{1}, c{2}, n);
%plot(sample(:,1), sample(:,2), '.');


