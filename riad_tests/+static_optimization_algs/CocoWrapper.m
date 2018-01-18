classdef CocoWrapper < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fHandle;
        fname;
        allEvals;
        invertSign;
        savePerfFunDir = [];
    end
    
    methods
        function [obj] = CocoWrapper(fHandle, fname, invertSign, savePerfFunDir)
            obj.fHandle = fHandle;
            obj.fname = fname;
            obj.invertSign = invertSign;
            obj.allEvals = [];
            if(exist('savePerfFunDir', 'var'))
                obj.savePerfFunDir = savePerfFunDir;
            end
        end
        function [vals] = eval(obj, x)
            %             DIM = size(x, 2);
            %             vals = (-obj.fHandle(x') / abs(obj.fHandle(zeros(DIM, 1))))';
            if(obj.invertSign)
                vals = -obj.fHandle(x')';
            else
                vals = obj.fHandle(x')';
            end
            obj.allEvals = [obj.allEvals; vals];
            if(~isempty(obj.savePerfFunDir) && obj.savePerfFunDir)
                currEvals = obj.allEvals;
                save(fullfile(obj.savePerfFunDir, ['currPerf' obj.getSignature]), 'currEvals');
            end
        end
        
        function [sign] = getSignature(obj)
            sign = obj.fname;
        end
        
        function plot(obj)
            %only valid if dim == 2
            limitInf = [-5, -5];
            limitSup = [5, 5];
            nbPoint = 99;
            line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
            line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
            [x1, x2] = meshgrid(line1, line2);
            y = obj.eval([x1(:) x2(:)]);
            contour(line1, line2, reshape(y, length(line2), length(line1)));
        end
    end
    
end

