function [xMin, yMin] = globalMinSearch(fcn, lb, ub)
    x0 = rand(1,size(lb,2)) .* (ub - lb) + lb;
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective',fcn,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    gs = GlobalSearch('Display', 'off'); %, 'MaxTime', 60);
    try
%         warning('off','MATLAB:singularMatrix');
        [xMin, yMin] = run(gs,problem);
%         warning('on','MATLAB:singularMatrix');
%         surfFcn2D(fcn,lb,ub,100);
%         hold on;
%         plot3(xMin(1), xMin(2), yMin,'ro');
%         xlabel('log(Hyper1)');
%         ylabel('log(Hyper2)');
%         zlabel('marginal likelihood');  
%         hold off;
        
        if any(xMin == lb) || any(xMin == ub)
            warning('off','backtrace');
            warning('Hyper parameters: minimization succeeded at boundary.');
            warning('on','backtrace');
        end
    catch me
        xMin = [];
        yMin = [];
        warning('off','backtrace');
        warning('Hyper parameters: minimization of marginal likelihood failed.');
        warning('on','backtrace');
    end 
end

