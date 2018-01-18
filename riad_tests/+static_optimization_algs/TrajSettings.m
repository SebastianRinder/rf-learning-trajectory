classdef TrajSettings < Common.IASObject
   properties (SetObservable, AbortSet)
       mean = 0;
       var = 1;
   end
   
   methods
       function [obj] = TrajSettings(mean, var)
           obj = obj@Common.IASObject();
           if(nargin > 0)
               obj.mean = mean;
           end
           if(nargin > 1)
               obj.var = var;
           end
           
           obj.linkProperty('mean');
           obj.linkProperty('var');
       end
   end
end
