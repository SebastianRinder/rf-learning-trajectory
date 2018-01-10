function hyper = optimizeHyper(knownX, y,D, opts)
%     hyperLb(1:2) = opts.hyper(1,1:2) - 3;
%     hyperUb(1:2) = opts.hyper(1,1:2) + 3;
%     
%     hyperLb(opts.hyper == 0) = -10;
%     hyperUb(opts.hyper == 0) = 10;

    hyperLb = [-10,-10];
    hyperUb = -hyperLb;

    negHyperFcn = @(X) log(-findHyper(X));
    [hyper, yMin] = globalMinSearch(negHyperFcn, hyperLb, hyperUb, opts.useGADSToolbox, true);

%     hyper(hyper < hyperLb + 0.5) = 0;
%     hyper(hyper > hyperUb - 0.5) = 0;
%     
%     if any(hyper == 0)
%         warning('off','backtrace');
%         warning('Hyper parameters: minimization succeeded near boundary.');
%         warning('on','backtrace');
%         yMin = negHyperFcn(hyper);
%     end    
%     hyper(1,hyper == 0) = opts.hyper(1,hyper == 0); %if a parameter is found near a boundary, take the last one

    margin = 0.5;
    
    hyper(hyper < hyperLb + margin) = 0;
    hyper(hyper > hyperUb - margin) = 0;

    if any(hyper == 0)
        warning('off','backtrace');
        warning('Hyper parameters: minimization succeeded near boundary.'); %maybe double like4
        warning('on','backtrace');
    end
    
    if opts.hyperPlot
        selectFigure('hyperPlot');
        clf;
        surfFcn2D(negHyperFcn,hyperLb,hyperUb,100);
        hold on;
        plot3(hyper(1), hyper(2), yMin,'ro','LineWidth',3);
        xlabel('hyper1: log(sigmaf)');
        ylabel('hyper2: log(sigmal)');
        zlabel('log(-marginal likelihood)');        

        title(['sigmaf: ',num2str(exp(hyper(1))),', sigmal: ',num2str(exp(hyper(2)))]);
        hold off;
        pause(0.1);
    end    

    function logLi = findHyper(hyper)
        %
%         opts.sigmaNoiseSquared = exp(hyper(1));
%         hyper(1) = 1;
        %
        K = opts.scaleKernel(D, hyper);
        L = getLowerCholesky(K, y, true, opts.sigmaNoiseSquared);

        if ~isempty(L)
            like1 = 2*sum(log(diag(L)));    
            like2 = y' * (L'\(L \ y));
            like3 = size(knownX,1) * log(2*pi);
            like4 = hyper(2).^2 ./ 200 - log(10*sqrt(2*pi)); %Lizotte 2008
            logLi = 0.5 * (-like1 - like2 - like3 + like4);                %Bishop 2006 (6.69)
        else
            logLi = -Inf;
        end
    end
end