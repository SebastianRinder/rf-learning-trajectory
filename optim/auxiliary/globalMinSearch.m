%% author: Sebastian Rinder

function [xMin, fMin] = globalMinSearch(fcn, lb, ub, useGADSToolbox, hyperSearch)
    if useGADSToolbox
        x0 = randBound(lb,ub,1);
        k = 0;
        while fcn(x0) == Inf && k < 1000
            x0 = randBound(lb,ub,1);
            k = k + 1;
        end
        if k >= 1000, error('no initial Point found'); end
        opts = optimoptions(@fmincon,'Algorithm','interior-point');
        problem = createOptimProblem('fmincon','objective',fcn,'x0',x0,'lb',lb,'ub',ub,'options',opts);
        gs = GlobalSearch('Display', 'off'); %, 'MaxTime', 60);

        [xMin, fMin] = run(gs,problem);
    else
        xRand = randBound(lb,ub,10000);
        if hyperSearch            
            for i=1:10000
                fRand(i,1) = fcn(xRand(i,:));
            end
        else
            fRand = fcn(xRand);
        end
        
        [~, rows] = sort(fRand, 'ascend');
        x0 = xRand(rows(1:10), :);
        
        options = optimoptions('fmincon','Display','off');
        
        warning('off','backtrace');
        warning off MATLAB:singularMatrix
        warning off MATLAB:nearlySingularMatrix
        for i = 1:10
            [x(i,:), f(i)] = fmincon(fcn,x0(i,:),[],[],[],[],lb,ub,[],options);
        end
        warning on MATLAB:singularMatrix
        warning on MATLAB:nearlySingularMatrix
        warning('on','backtrace');

        [fMin, idx] = min(f);
        xMin = x(idx,:);
    end
end
