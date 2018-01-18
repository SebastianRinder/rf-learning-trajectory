function [L, alpha] = getLowerCholesky(K, knownY, calcHyper, noiseVariance)    
    alpha = [];    
    
    p = 1;
    doublings = 0;
    while p > 0 && doublings < 100
        Knoise = K + noiseVariance.*eye(size(K));
        [L,p] = chol(Knoise, 'lower');
        if p > 0
            L = [];
            if ~calcHyper
                doublings = doublings + 1;
                noiseVariance = noiseVariance * 2;
            else
                break;
            end
        end        
    end
%     global nv;
%     nv = noiseVariance;

    if nargout > 1 && ~isempty(L)
        alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
    end
end