function D = trajectoryCovariance(Xm, Xn, trajectories, opts)
    m = size(Xm,1);
    n = size(Xn,1);
    
    D = zeros(m,n);    
    
    prior = false;
    if m == n && all(all(Xm == Xn))
        prior = true;
    end
    
    for i = 1:m
        start = 1;
        if prior, start = i; end
        for j = start:n
            if all(Xm(i,:) == Xn(j,:))
                D(i,j) = 0;
            else
%                 traji = findTraj(Xm(i,:), trajectories);
%                 trajj = findTraj(Xn(j,:), trajectories);
                if prior
                    D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), trajectories(i,:), trajectories(j,:), opts);
                else
                    D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), [], trajectories(j,:), opts);
                end
            end
        end
    end
    
    if prior
        D = D + D';
    end
    
end


function traj = findTraj(X, trajectory)
    idx = find(all(X == trajectory.policy, 2));
    if isempty(idx)
        traj = [];
    else
        traj = trajectory.data(idx,:);
    end
end