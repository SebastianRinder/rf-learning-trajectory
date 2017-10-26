function [mean, variance] = gaussianProcessModel(testX, knownX, knownY, preComputed, hyper)
    Ks = sqExpCovariance(testX, knownX, hyper);
    
    if nargout == 1
        alpha = preComputed;
        mean = Ks * alpha;      %Rasmussen, Williams 2005
    else
        L = preComputed;    
        Kss = 1;
        
        Lk = L \ Ks';       %Python example, Nando de Freitas 2013
        mean = Lk' * (L \ knownY);
        variance = (Kss - sum(Lk.^2))';
    end
end