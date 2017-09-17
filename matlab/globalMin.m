function [x,f] = globalMin(fun,x0,lb,ub)
    x0 = randTheta(lb,ub);
    %opts = optimoptions(@fmincon,'Algorithm','interior-point');
    %problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',lb,'ub',ub,'options',opts);
    %gs = GlobalSearch('Display', 'iter'); %, 'MaxTime', 60);
    %[x,f] = run(gs,problem);
    
    %[x,f] = fmincon(fun,x0,[],[],[],[],lb,ub);
    
    opts = optimoptions(@patternsearch,'MeshTolerance',1e-6);
    [x,f] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],opts);
end

