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

            ti = find(all(Xm(i,:) == trajectory.policy, 2));
            if isempty(ti)
                traji = [];
            else
                if any(size(ti) ~= [1 1])
                    keyboard;
                end
                traji = trajectory.data(ti,:);
            end

            tj = find(all(Xn(j,:) == trajectory.policy, 2));
            if isempty(tj)
                trajj = [];
            else
                if any(size(tj) ~= [1 1])
                    keyboard;
                end
                trajj = trajectory.data(tj,:);
            end

            D(i,j) = trajectoryDistance(Xm(i,:), Xn(j,:), traji, trajj, customOpts);
        end
    end
    
    K = scale(D, hyper);
end

function K = scale(D, hyper)
%     m = max(max(D));
%     if m == 0
%         K = hyper.f .* ones(size(D));
%     else
%         hyper.l = 10/m;
%         K = hyper.f .* exp(-hyper.l .* D);
%     end
    K = hyper.f .* exp(-hyper.l .* D);
end