function K = trajectoryCovariance(Xm, Xn, theta, customOpts)
    if ~isstruct(theta)
        hyper.l = theta(1);
        hyper.f = theta(2);
    else
        hyper = theta;
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
    K = scale(D, hyper);
end

function traj = findTraj(X, trajectory)
    idx = find(all(X == trajectory.policy, 2));
    if isempty(idx)
        traj = [];
    else
        traj = trajectory.data(idx,:);
    end
end

function K = scale(D, hyper)
    hyper.f = 1;

    if max(max(D)) == 0
        K = hyper.f .* ones(size(D));
    else
        m = min(min(D(D~=0)));
        if m >= 10
            hyper.l = 0.1/m;
        else
            hyper.l = 1;
        end
        K = hyper.f .* exp(-hyper.l .* D);
    end
end