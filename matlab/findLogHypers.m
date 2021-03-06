function logLi = findLogHypers(logHyper, knownX, knownY, opts)
    opts.hyper = exp(logHyper);
    
    L = getLowerCholesky(opts.D, knownY, opts, true);
    
    try
        if ~isempty(L)
            like1 = 2*sum(log(diag(L)));    
            like2 = knownY' * (L'\(L \ knownY));
            like3 = size(knownX,1) * log(2*pi);
            logLi = -0.5 * (like1 + like2 + like3); %Bishop 2006 (6.69)
        else
            logLi = -Inf;
        end
    catch me
        disp('fail');
    end
end