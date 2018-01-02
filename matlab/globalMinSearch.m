function [xMin, yMin] = globalMinSearch(fcn, lb, ub,toHyper,~)
    x0 = randPolicy(lb,ub,1);
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective',fcn,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    gs = GlobalSearch('Display', 'off'); %, 'MaxTime', 60);
    [xMin, yMin] = run(gs,problem);
end
