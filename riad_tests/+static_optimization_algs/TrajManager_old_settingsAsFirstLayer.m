classdef TrajManager < Data.DataManager
    methods
        function obj = TrajManager(sizeContext, sizeParams, sizeState)
            obj = obj@Data.DataManager('meta');

            obj.addDataEntry('mean', 1);
            obj.addDataEntry('var', 1);
            
            obj.addDataAlias('metaParams', {'mean', 'var'});

            
            trajManager = Data.DataManager('trajs');

            
            trajManager.addDataEntry('context', sizeContext);
            trajManager.addDataEntry('params', sizeParams);
            
            trajManager.addDataAlias('trajConf', {'context', 'params'});

            subManager = Data.DataManager('states');
            subManager.addDataEntry('state', sizeState);
            trajManager.setSubDataManager(subManager);
            
            obj.setSubDataManager(trajManager);
            
            obj.finalizeDataManager();
            
        end
    end
end
