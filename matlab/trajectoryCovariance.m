function Kmn = trajectoryCovariance(Xm, Xn, theta, customOpts)
    if ~isstruct(theta)
        hyper.l = max(theta(1), 1e-6);
        hyper.f = max(theta(2), 1e-6);
    else
%         hyper = theta;
    end

    m = size(Xm,1);
    n = size(Xn,1);
    
    D = zeros(m,n);    
    
    prior = false;
    if m == n && all(all(Xm == Xn))
        prior = true;
    end
    
    for i = 1:m
        start = 1;
        if prior
            start = i;
        end
        for j = start:n
            if all(Xm(i,:) == Xn(j,:))
                D(i,j) = 0;
            else
                traji = findTraj(Xm(i,:), customOpts.trajectory);
                trajj = findTraj(Xn(j,:), customOpts.trajectory);
                D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), traji, trajj, customOpts);
            end
        end
    end
    
    if prior
        D = D + D';
    end    
    Kmn = scale(D, hyper);
end

function traj = findTraj(X, trajectory)
    idx = find(all(X == trajectory.policy, 2));
    if isempty(idx)
        traj = [];
    else
        traj = trajectory.data(idx,:);
    end
end

function Kmn = scale(D, hyper)
    if max(max(D)) == 0
        Kmn = hyper.f .* ones(size(D));
    else
        m = 16 / max(max(D));
        d = m * D;
        
        Kmn = d ./ (hyper.l^2);
        Kmn = (hyper.f^2) * exp(-0.5*Kmn);
    end
end