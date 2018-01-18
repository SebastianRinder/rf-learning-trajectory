classdef TestMapping < Functions.Mapping
   methods
       function obj = TestMapping(dataManager)
           obj = obj@Functions.Mapping(dataManager, 'sumparam', 'params', 'testMap');       
           
           obj.addMappingFunction('Sum');
           obj.addMappingFunction('NR', 'notRegistred');
       end
       
       function res = Sum(obj, numElements, params)
           res = sum(params, 2);       
       end
       
       function res = NR(obj, numElements, params)
           res = mod(sum(params, 2), 4);       
       end
       
   end
end