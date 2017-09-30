function K = trajectoryCovarianceOLD(Xm, Xn, trajectory, hyper, actionSelectionFun)
    if isempty(Xm) %compute K from prior
        K = zeros(size(Xn,1));
        for i = 1:size(Xn,1)
            for j = i:size(Xn,1)
                if i==j
                    K(i,j) = 0.5; %adding the transposed matrix later sums the diagonal elements to 1
                else
                    K(i,j) = kernelTrajectory(Xn(i,:), Xn(j,:), trajectory(i,:), trajectory(j,:), hyper, actionSelectionFun);
                end
            end
        end
        K = K + K'; %K is symmetric because of symmetric kernel
    else %compute K from prior and posterior
        K = zeros(1,size(Xn,1));
        for j = 1:size(Xn,1)
            K(1,j) = kernelTrajectory(Xm, Xn(j,:), [], trajectory(j,:), hyper, actionSelectionFun);
        end
    end
end