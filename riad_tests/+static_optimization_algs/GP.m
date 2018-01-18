classdef GP
    methods (Static)
        function addGPStuffPath(gpstuffroot)
            addpath([gpstuffroot 'diag'])
            addpath([gpstuffroot 'dist'])
            addpath([gpstuffroot 'gp'])
            addpath([gpstuffroot 'mc'])
            addpath([gpstuffroot 'misc'])
            addpath([gpstuffroot 'optim'])
            addpath([gpstuffroot 'tests'])
        end
        
        function plotGP(gp, xl, yl, center, width, mu, rotMatrix, featureFunction)
            dim = length(center);
            if(dim == 2)
                nbPoint = 100;
                limitInf = center - width;
                limitSup = center + width;
                line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
                line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
                [x1, x2] = meshgrid(line1, line2);
                x = [x1(:) x2(:)];
                if(exist('mu', 'var') && ~isempty(mu))
                    x = bsxfun(@minus, x, mu);
                end
                if(exist('rotMatrix', 'var') && ~isempty(rotMatrix))
                    x = x * rotMatrix;
                end
                if(exist('featureFunction', 'var') && ~isempty(featureFunction))
                    x = featureFunction(x);
                end
                
                [y, ~] = gpmc_pred(gp, xl, yl, x);
                %                 h = surf(line1, line2, reshape(y, length(line2), length(line1)));
                %                 set(h,'LineStyle','none');
                contour(line1, line2, reshape(y, length(line2), length(line1)));
            else
                warning('plot not implemented for this dimension');
            end
        end
        
        function plotGPML(hyp_param, infGaussLik, meanfunc, covfunc, likfunc, xl, yl, center, width, mu, rotMatrix, featureFunction)
            dim = length(center);
            if(dim == 2)
                nbPoint = 100;
                limitInf = center - width;
                limitSup = center + width;
                line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
                line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
                [x1, x2] = meshgrid(line1, line2);
                x = [x1(:) x2(:)];
                if(exist('mu', 'var') && ~isempty(mu))
                    x = bsxfun(@minus, x, mu);
                end
                if(exist('rotMatrix', 'var') && ~isempty(rotMatrix))
                    x = x * rotMatrix;
                end
                if(exist('featureFunction', 'var') && ~isempty(featureFunction))
                    x = featureFunction(x);
                end
                
                y = gp(hyp_param, infGaussLik, meanfunc, covfunc, likfunc, xl, yl, x);
                %                 h = surf(line1, line2, reshape(y, length(line2), length(line1)));
                %                 set(h,'LineStyle','none');
                contour(line1, line2, reshape(y, length(line2), length(line1)));
            else
                warning('plot not implemented for this dimension');
            end
        end
        
        function plotGPRand(gp, xl, yl, center, width, mu, rotMatrix, featureFunction)
            dim = length(center);
            if(dim == 2)
                nbPoint = 100;
                limitInf = center - width;
                limitSup = center + width;
                line1 = limitInf(1):(limitSup(1)-limitInf(1))/nbPoint:limitSup(1);
                line2 = limitInf(2):(limitSup(2)-limitInf(2))/nbPoint:limitSup(2);
                [x1, x2] = meshgrid(line1, line2);
                x = [x1(:) x2(:)];
                if(exist('mu', 'var') && ~isempty(mu))
                    x = bsxfun(@minus, x, mu);
                end
                if(exist('rotMatrix', 'var') && ~isempty(rotMatrix))
                    x = x * rotMatrix;
                end
                if(exist('featureFunction', 'var') && ~isempty(featureFunction))
                    x = featureFunction(x);
                end
                [y, ~] = gp_rnd(gp, xl, yl, x);
                %                 h = surf(line1, line2, reshape(y, length(line2), length(line1)));
                %                 set(h,'LineStyle','none');
                contour(line1, line2, reshape(y, length(line2), length(line1)));
            else
                warning('plot not implemented for this dimension');
            end
        end
        
        function [m, v] = predFromCov(gp, xt, x, Cy, invCholC)
            kxx = gp_trcov(gp, xt);
            kx = gp_cov(gp, xt, x)';
            m = kx' * Cy;
            rootKC = invCholC' * kx;
            v = kxx - (rootKC' * rootKC);
        end
        
        function val = randFromCov(gp, xt, x, Cy, invCholC)
            [m, v] = static_optimization_algs.GP.predFromCov(gp, xt, x, Cy, invCholC);
            val = m + randn * sqrt(v);
%             val = m;
        end

        function val = EI(gp, xt, x, Cy, invCholC, target)
            [m, v] = static_optimization_algs.GP.predFromCov(gp, xt, x, Cy, invCholC);
            sv = sqrt(v);
            val = 0;
            if(sv > 1e-16)
                z = (m - target) / sv;
                val = sv * (z * normcdf(z) + normpdf(z));            
            end                
        end        

        function val = EILib(gp, xt, x, y, target)
            [m, v] = gp_pred(gp, x, y, xt);
            sv = sqrt(v);
            val = 0;
            if(sv > 1e-16)
                z = (m - target) / sv;
                val = sv^2 * (z * normcdf(z) + normpdf(z));
            end
        end
        
        
        function val = ThompsonSampling(gp, xt, x, y)
            persistent xprev;
            persistent yprev;
            if(isempty(xt))
                xprev = [];
                yprev = [];
                return;
            end            
            val = gp_rnd(gp, [x; xprev], [y; yprev], xt);
            xprev = [xprev; xt];
            yprev = [yprev; val];
        end

        function val = UCB(gp, xt, x, Cy, invCholC, beta)
            [m, v] = static_optimization_algs.GP.predFromCov(gp, xt, x, Cy, invCholC);
            val = m + beta * sqrt(v);
        end
        
        function xstar = localArgmaxSample(gp, xl, yl, normal, nbSampleGauss, mu, rotMatrix, featureFunction)
            samples = normal.getSamples(nbSampleGauss);
            transSamples = samples;
            if(exist('mu', 'var') && ~isempty(mu))
                transSamples = bsxfun(@minus, transSamples, mu);
            end
            if(exist('rotMatrix', 'var') && ~isempty(rotMatrix))
                transSamples = transSamples * rotMatrix;
            end
            if(exist('featureFunction', 'var') && ~isempty(featureFunction))
                transSamples = featureFunction(transSamples);
            end
            [ysamples, ~] = gp_rnd(gp, xl, yl, transSamples);
            [~, ixstar] = max(ysamples);
            xstar = samples(ixstar, :)';
        end
        
        function vals = gpRandTrans(gp, xl, yl, samples, mu, rotMatrix, featureFunction)
            if(exist('mu', 'var') && ~isempty(mu))
                samples = bsxfun(@minus, samples, mu);
            end
            if(exist('rotMatrix', 'var') && ~isempty(rotMatrix))
                samples = samples * rotMatrix;
            end
            if(exist('featureFunction', 'var') && ~isempty(featureFunction))
                samples = featureFunction(samples);
            end
            vals = gp_rnd(gp, xl, yl, samples);
        end
      
        function vals = gpPredTrans(gp, xl, yl, samples, mu, rotMatrix, featureFunction)
            if(exist('mu', 'var') && ~isempty(mu))
                samples = bsxfun(@minus, samples, mu);
            end
            if(exist('rotMatrix', 'var') && ~isempty(rotMatrix))
                samples = samples * rotMatrix;
            end
            if(exist('featureFunction', 'var') && ~isempty(featureFunction))
                samples = featureFunction(samples);
            end
            vals = gp_pred(gp, xl, yl, samples);
        end

        
    end
end

