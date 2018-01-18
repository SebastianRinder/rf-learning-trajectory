classdef TrajManipulator < Data.DataManipulator
    
    properties (Access = protected)
        trajSettings;
    end
    
    methods
        function obj = TrajManipulator(dataManager)
            obj = obj@Data.DataManipulator(dataManager);
            
            obj.trajSettings = Common.Settings();
            
            obj.addDataManipulationFunction('perturbParams', {'params'}, {'params'});
            
            obj.addDataManipulationFunction('updateContext', {'params', 'context'}, {'context'});
            
            obj.addDataManipulationFunction('filterStates', {'context', 'state'}, {'state'}, Data.DataFunctionType.PER_EPISODE);
            
            obj.addDataManipulationFunction('sumStates', {'state'}, {});

        end
        
        function [params] = perturbParams(obj, params)
            params = params + obj.trajSettings.getProperty('mean') + obj.trajSettings.getProperty('var') * randn(size(params));
        end
        
        function [context] = updateContext(obj, params, context)
            context = round((context + params) .* rand(size(context)));
        end
        
        function [state] = filterStates(obj, context, state)
            state = state .* abs(mod(sum(context), 5) - mod(state, 5));
        end
        
        function [res] = sumStates(obj, state)
           res = sum(state, 2); 
        end 
    end
    
end
