clear all;
clc;

dataManager = Data.DataManager('dataset');

dataManager.addDataEntry('x', 5);
dataManager.addDataEntry('y', 2);
dataManager.addDataEntry('alpha', 1);
myData = dataManager.getDataObject(1000);

myData.setDataEntry('x', rand(myData.getNumElements('x'), myData.getNumDimensions('x')));

alpha = (1:myData.getNumElements('x'))';
myData.setDataEntry('alpha', alpha);

weights = rand(myData.getNumDimensions('y'), myData.getNumDimensions('x'));
bias = rand(myData.getNumDimensions('y'),1);
temp = rand(myData.getNumDimensions('y'), myData.getNumDimensions('y'));
chola = chol(temp' * temp);

gaussianDistrib = Distributions.Gaussian.GaussianLinearInFeatures(dataManager, 'y', 'x', 'gaussianDistrib');
gaussianDistrib.initObject();

gaussianDistrib.setWeightsAndBias(weights, bias);
gaussianDistrib.setSigma(chola);

gaussianDistrib.callDataFunction('sampleFromDistribution', myData);

learner = Learner.SupervisedLearner.LinearGaussianMLLearner(dataManager, gaussianDistrib);

sum(sum(abs(weights - learner.functionApproximator.weights)))

learner.updateModel(myData);

sum(sum(abs(weights - learner.functionApproximator.weights)))

learner.setWeightName('alpha');

learner.updateModel(myData);

sum(sum(abs(weights - learner.functionApproximator.weights)))


