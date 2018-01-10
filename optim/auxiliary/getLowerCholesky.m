function [L, alpha] = getLowerCholesky(K, knownY, calcHyper, sigmaNoiseSquared)    
    alpha = [];    
    
    p = 1;
    doublings = 0;
    while p > 0 && doublings < 100
        Knoise = K + sigmaNoiseSquared.*eye(size(K));
        [L,p] = chol(Knoise, 'lower');
        if p > 0
            L = [];
            if ~calcHyper
                doublings = doublings + 1;
                sigmaNoiseSquared = sigmaNoiseSquared * 2;
            else
                break;
            end
        end        
    end
    
%     if doublings > 0
%         disp(sigmaNoiseSquared);
%     end
    
    if nargout > 1 && ~isempty(L)
        alpha = L'\(L \ knownY);        %Rasmussen, Williams 2005
    end
end