classdef StateFeatures < FeatureGenerators.FeatureGenerator
    
    methods
        function [obj] = StateFeatures(dataManager, featuresIn)
            numFeatures = dataManager.getNumDimensions('state');
            numFeatures = numFeatures * numFeatures + 1;
            obj = obj@FeatureGenerators.FeatureGenerator(dataManager, featuresIn, 'Squared', ':', numFeatures);            
        end
        
        function [features] = getFeaturesInternal(obj, ~, inputMatrix)
            features = [ones(size(inputMatrix, 1), 1) (inputMatrix) (inputMatrix .* inputMatrix)];
        end
        
    end
    
end
