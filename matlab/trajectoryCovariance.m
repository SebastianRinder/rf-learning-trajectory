function K = trajectoryCovariance(Xm, Xn, theta, traj, actionSelectionFcn)
    if ~isstruct(theta)
        hyper.l = theta(1);
        hyper.f = theta(2);
    else
        hyper = theta;
    end
    
    calcPosterior = true;
    
    m = size(Xm,1);
    n = size(Xn,1);
    if m == n || n == 1 && all(all(Xm == Xn))
        calcPosterior = false;
    end
    K = zeros(m,n);
    
    if calcPosterior
        for i = 1:m
            for j = 1:n
                K(i,j) = kernelTrajectory(Xm(i,:), Xn(j,:), [], traj(j,:), hyper, actionSelectionFcn);
            end
        end
    else
        for i = 1:m
            for j = i:n
                if i==j
                    K(i,j) = hyper.f/2; %adding the transposed matrix later sums the diagonal elements to hyper.f
                else
                    K(i,j) = kernelTrajectory(Xm(i,:), Xn(j,:), traj(i,:), traj(j,:), hyper, actionSelectionFcn);
                end
            end
        end
        K = K + K'; %K is symmetric because of symmetric kernel
    end
end