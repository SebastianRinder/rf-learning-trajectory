classdef cartPoleClass < handle
    properties
        opts;
    end
    
    methods
        function obj = cartPoleClass()
            obj.opts = environmentSettings('cartPole', 'none');
            obj.opts.trajectoriesPerPolicy = 5;
        end
        
        function vals = eval(obj, samples)
            vals = zeros(size(samples,1),1);
            for i=1:size(samples,1)
                vals(i,1) = objectiveFcn(samples(i,:), obj.opts);
            end
        end
    end
    
end