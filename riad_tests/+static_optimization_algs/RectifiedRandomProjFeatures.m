classdef RectifiedRandomProjFeatures < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        W;
        b;
    end
    
    methods
        function [obj] = RectifiedRandomProjFeatures(dim, nbFeatures, sigParams)
            obj.W = sigParams * randn(dim, nbFeatures);
            obj.b = sigParams * randn(1, nbFeatures);
        end
        
        function fx = getFeatures(obj, x)
            fx = [ones(size(x, 1), 1) log( 1 + exp( - abs( bsxfun(@plus, x * obj.W, obj.b) ) ) )];
        end        
    end    
end
