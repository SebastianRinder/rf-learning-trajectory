function D = trajectoryDistanceLocal(Xm, Xn, trajectories, isThompsonsSample, opts)
    m = size(Xm,1);
    n = size(Xn,1);
    
    D = zeros(m,n);    
    
    states = [];
    for i = 1:size(trajectories,1) %only 1 sample per trajectory
        states = [states; trajectories{i,1}.state];
    end
    if size(states,1) > 500
        states = states(randsample(size(states,1),500),:);
    end
    
    if m == n && all(all(Xm == Xn)) %symmetric = true
        for i = 1:m
            for j = i:n
                if i==j
                    D(i,j) = 0;
                else
                    D(i,j) = stateDistance(Xm(i,:), Xn(j,:), states, opts);
                end
            end
        end

        D = D + D';
    else
        for i = 1:m
            for j = 1:n
                D(i,j) = stateDistance(Xm(i,:), Xn(j,:), states, opts);
            end
        end
    end
end

function D = stateDistance(xi, xj, states, opts)
    if isequal(opts.actionSpace, 'continuous')
        [~,~,mu1] = opts.actionSelectionFcn(xi, states, [], opts.actionMisc);
        [~,~,mu2] = opts.actionSelectionFcn(xj, states, [], opts.actionMisc);
        muDiff = mu1 - mu2;
        D = (muDiff' * muDiff) ./ size(states,1);
    else
        [~,~,prob1] = opts.actionSelectionFcn(xi, states, [], opts.actionMisc);
        [~,~,prob2] = opts.actionSelectionFcn(xj, states, [], opts.actionMisc);
        probDiff = sum(abs(prob1 - prob2),2);
%         D = (probDiff' * probDiff) ./ size(states,1); %take mean instead
        D = mean(probDiff);
    end
end