function [xMin, yMin] = globalMinSearch(fcn, lb, ub)
    x0 = rand(1,size(lb,2)) .* (ub - lb) + lb;
    k = 0;
    while fcn(x0) == Inf && k < 1000
        x0 = rand(1,size(lb,2)) .* (ub - lb) + lb;
        k = k + 1;
    end
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective',fcn,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    gs = GlobalSearch('Display', 'off'); %, 'MaxTime', 60);
   try
%         warning('off','MATLAB:singularMatrix');
        [xMin, yMin] = run(gs,problem);
%         warning('on','MATLAB:singularMatrix');
        surfFcn2D(fcn,lb,ub,100);
        hold on;
        plot3(xMin(1), xMin(2), yMin,'ro','LineWidth',3);
        xlabel('log(hyper1)');
        ylabel('log(hyper2)');
        zlabel('log(-marginal likelihood)');  
                
        xMin(xMin < lb + 1) = 0;
        xMin(xMin > ub - 1) = 0;
        
        title(['hyper1: ',num2str(exp(xMin(1))),', hyper2: ',num2str(exp(xMin(2)))]);
        hold off;
        pause(0.1);
        
        if any(xMin == 0)
            warning('off','backtrace');
            warning('Hyper parameters: minimization succeeded near boundary.');
            warning('on','backtrace');
            yMin = [];
        end
    catch me
        keyboard;
%         xMin = [];
%         yMin = [];
%         warning('off','backtrace');
%         warning('Hyper parameters: minimization of marginal likelihood failed.');
%         warning('on','backtrace');
    end 
end

