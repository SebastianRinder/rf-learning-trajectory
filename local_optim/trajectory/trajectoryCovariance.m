function D = trajectoryCovariance(Xm, Xn, trajectories, isCovMat, opts)
    m = size(Xm,1);
    n = size(Xn,1);
    
    D = zeros(m,n);    
    
    if m == n && all(all(Xm == Xn)) %symmetric = true
        if ~isCovMat
            for i = 1:m
                for j = i:n
                    if all(Xm(i,:) == Xn(j,:))
                        D(i,j) = 0;
                    else
                        D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), trajectories(i,:), trajectories(j,:), opts);
                    end
                end
            end
        else
            for i = 1:m
                for j = i:n
                    if all(Xm(i,:) == Xn(j,:))
                        D(i,j) = 0;
                    else
                        D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), trajectories, [], opts);
                    end
                end
            end
        end
        D = D + D';
        
    else
        for i = 1:m
            for j = 1:n
                if all(Xm(i,:) == Xn(j,:))
                    D(i,j) = 0;
                else
                    D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), [], trajectories(j,:), opts);
                end
            end
        end
    end
end