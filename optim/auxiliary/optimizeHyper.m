function hyper = optimizeHyper(knownX, y,D, opts)
    lb = [-10,opts.hyper(2)-10];
    ub = [10,opts.hyper(2)+10];
%     lb = opts.hyperLb;
%     ub = opts.hyperUb;
    center = (ub - lb)./2 + lb;
    
    negHyperFcn = @(X) log(-findHyper(X));
    [hyper, yMin] = globalMinSearch(negHyperFcn, lb, ub, opts.useGADSToolbox, true);
    if isempty(hyper)
        hyper = opts.hyper;
    else
        if opts.hyperPlot && size(lb,2) == 2
            selectFigure('hyperPlot');
            clf;
            plotFcn(negHyperFcn,lb,ub,100,2);
            hold on;
            plot3(hyper(1), hyper(2), yMin,'ro','LineWidth',3);
            xlabel('hyper1: log(sigmaf)');
            ylabel('hyper2: log(sigmal)');
            zlabel('log(-marginal likelihood)');        

            title(['sigmaf: ',num2str(exp(hyper(1))),', sigmal: ',num2str(exp(hyper(2)))]);
            hold off;
            pause(0.1);
        elseif opts.hyperPlot && size(lb,2) == 1
            selectFigure('hyperPlot');
            clf;
            plotFcn(negHyperFcn,lb,ub,100,1);
            hold on;
            plot(hyper, yMin,'ro','LineWidth',3);
            xlabel('hyper2: log(sigmal)');
            ylabel('log(-marginal likelihood)');        

            title(['sigmal: ',num2str(exp(hyper))]);
            hold off;
            pause(0.1);
        end

        margin = 0.5;        
        hyper(hyper < lb + margin) = 0;
        hyper(hyper > ub - margin) = 0;

        if any(hyper == 0)
            warning('off','backtrace');
            warning('Hyper parameters: minimization succeeded near boundary.');
            warning('on','backtrace');
            %hyper(hyper == 0) = opts.hyper(hyper==0);
            hyper(hyper == 0) = center(hyper==0);
        end 
    end

    function logLi = findHyper(hyper)
        if size(hyper,2) == 1
            K = opts.scaleKernel(D, [0,hyper]);
        else
            K = opts.scaleKernel(D, hyper);
        end
        L = getLowerCholesky(K, y, true, opts.noiseVariance);

        if ~isempty(L)
            like1 = 2*sum(log(diag(L)));    
            like2 = y' * (L'\(L \ y));
            like3 = size(knownX,1) * log(2*pi);
%             like4 = sum(-(hyper-center).^2 ./ (2.*(ub-lb).^2) - log((ub-lb)*sqrt(2*pi))); %Lizotte 2008 & move center
%             logLi = 0.5 * (-like1 - like2 - like3 + like4);                
            logLi = 0.5 * (-like1 - like2 - like3);
        else
            logLi = -Inf;
        end
    end
end