function K = trajectoryCovarianceCartPole(thetai, thetaj, hyper, trajectory)
    if isempty(thetaj) %compute K from prior
        K = zeros(size(thetai,1));
        for i = 1:size(thetai,1)
            for j = i:size(thetai,1)
                if i==j
                    K(i,j) = 0.5; %adding the transposed matrix later sums the diagonal elements to 1
                else
                    K(i,j) = kernelKnownTrajectoryCartPole(thetai(i,:), thetai(j,:), trajectory(i,:), trajectory(j,:), hyper);
                end
            end
        end
        K = K + K'; %K is symmetric because of symmetric kernel
    else %compute Ks from prior and posterior
        K = zeros(1,size(thetaj,1));
        for j = 1:size(thetaj,1)
            K(1,j) = kernelUnknownTrajectoryCartPole(thetai, thetaj(j,:), trajectory(j,:), hyper);
        end
    end
end

