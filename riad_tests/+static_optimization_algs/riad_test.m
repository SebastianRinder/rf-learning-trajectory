clear all;
close all;
clc;

% dataManager = Data.DataManager('traj');
% dataManager.addDataEntry('context', 5);
% dataManager.addDataEntry('params', 5);
% 
% subManager = Data.DataManager('states');
% subManager.addDataEntry('state', 2);
% dataManager.setSubDataManager(subManager);
% dataManager.addDataAlias('trajConf', {'context', 'params'});
% 
% dataManager.finalizeDataManager();

contextDim = 5;
paramsDim = 5;
stateDim = 2;

trajManager = TrajManager(contextDim, paramsDim, stateDim); %dim of: context, params, states 

squared = StateFeatures(trajManager, 'state');

trajManager.getDataStructure().states

nbTrajs = 4;
statesPerTraj = 6;

myData = trajManager.getDataObject([nbTrajs statesPerTraj]);

meanGauss = 10;
varGauss = 2;

myXpSettings = TrajSettings(meanGauss, varGauss);

%Filling the data with sequential numbers
tempDim = [myData.getNumElements('context') myData.getNumDimensions('context')];
tempSize = tempDim(1) * tempDim(2);
rawData = reshape(1:tempSize, tempDim(1), tempDim(2));
myData.setDataEntry('context', rawData);

for i=1:myData.getNumElements('context')
    tempDim = [myData.getNumElements('state') myData.getNumDimensions('state')];
    rawData = reshape((tempSize + 1):(tempSize + tempDim(1) * tempDim(2)), tempDim(1), tempDim(2));
    tempSize = tempSize + tempDim(1) * tempDim(2);
    myData.setDataEntry('state', rawData, i);
end



myData.getDataEntry('context')
myData.getDataEntry('context', 2)
myData.getDataEntry('state')
myData.getDataEntry('state', 3)
myData.getDataEntry('state', 3, 4)
myData.getDataEntry('trajConf')
cellArray = myData.getDataEntryCellArray({'context', 'state'}, :, 2)
cellArray{1}
cellArray{2}

%ok now I am handling the data pretty ok... 
%let's see these DataManipulators
disp('----------Manipulator--------')

myManip = TrajManipulator(trajManager);

disp('Traj befor calling perturbParams')
myData.getDataEntry('trajConf')

myManip.callDataFunction('perturbParams', myData);

disp('Traj after calling perturbParams')
myData.getDataEntry('trajConf')

myManip.callDataFunction('filterStates', myData)

% for i=1:myData.getNumElements('context')
%     myManip.callDataFunction('filterStates', myData, i);
% end

myData.getDataEntry('state')

myManip.callDataFunctionOutput('sumStates', myData)

myManip.addDataFunctionAlias('perturbSum', 'perturbParams');
myManip.addDataFunctionAlias('perturbSum', 'updateContext');
myManip.addDataFunctionAlias('perturbSum', 'filterStates');
myManip.addDataFunctionAlias('perturbSum', 'sumStates');

for i=1:0
    myManip.callDataFunction('updateContext', myData);
    myData.getDataEntry('context')
end

% okiiiii let's do Settings now.. done, see above

% Features, anyone?
size(myData.getDataEntry('stateSquaredTag'))
myData.getDataEntry('stateSquared')

myManip.callDataFunction('filterStates', myData);

myData.getDataEntry('stateSquared')
myData.setDataEntry('stateSquaredTag', zeros(myData.getNumElements('stateSquaredTag'), 1));

myData.getDataEntry('stateSquared')

% Mappings
myMapp = TestMapping(trajManager);
myData.getDataEntry('sumparam')
myMapp.callDataFunction('Sum', myData);
myData.getDataEntry('sumparam')
myMapp.callDataFunctionOutput('NR', myData)

% Distributions now

