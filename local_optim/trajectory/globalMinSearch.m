function [xMin, yMin] = globalMinSearch(fcn, lb, ub)
    x0 = rand(1,size(lb,2)) .* (ub - lb) + lb;
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective',fcn,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    gs = GlobalSearch('Display', 'off'); %, 'MaxTime', 60);
    try
        [xMin, yMin] = run(gs,problem);
    catch me
%         surfFcn2D(fcn,lb,ub,100);
        disp('Hyper: no min found');
        
        xMin = [];
        yMin = [];
    end
    surfFcn2D(fcn,lb,ub,100);
end