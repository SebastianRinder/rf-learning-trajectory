function K = trajectoryCovariance(Xm, Xn, theta, customOpts)

    trajectory = customOpts.trajectory;
    
    if ~isstruct(theta)
        hyper.l = theta(1);
        hyper.f = theta(2);
    else
        hyper = theta;
    end

    m = size(Xm,1);
    n = size(Xn,1);
    
    D = zeros(m,n);    
    for i = 1:m
        for j = 1:n
            
            if all(Xm(i,:) == Xn(j,:))
                D(i,j) = 0;
            else
                ti = find(all(Xm(i,:) == trajectory.policy, 2));
                if isempty(ti)
                    traji = [];
                else
                    traji = trajectory.data(ti,:);
                end
                
                traji = find(Xm(i,:), trajectory);
                trajj = find(Xn(j,:), trajectory);

                tj = find(all(Xn(j,:) == trajectory.policy, 2));
                if isempty(tj)
                    trajj = [];
                else
                    trajj = trajectory.data(tj,:);
                end
                D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), traji, trajj, customOpts);
            end            
        end
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

%85 iters, 11705s, best at 38, best -0.149

function K = scale(D, hyper)
%     m = max(max(D));
%     hyper.f = 1;
%     if m == 0
%         K = hyper.f .* ones(size(D));
%     else
%         hyper.l = 10/m;
%         K = hyper.f .* exp(-hyper.l .* D);
%     end
    K = hyper.f .* exp(-hyper.l .* D);
end