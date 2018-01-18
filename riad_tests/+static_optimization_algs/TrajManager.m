classdef TrajManager < Data.DataManager
    methods
        function obj = TrajManager(sizeContext, sizeParams, sizeState)
            obj = obj@Data.DataManager('trajs');
            
            obj.addDataEntry('context', sizeContext);
            obj.addDataEntry('params', sizeParams);
            obj.addDataEntry('sumparam', 1);
            
            obj.addDataAlias('trajConf', {'context', 'params'});

            subManager = Data.DataManager('states');
            subManager.addDataEntry('state', sizeState);
            obj.setSubDataManager(subManager);
            
            obj.finalizeDataManager();
            
        end
    end
end
