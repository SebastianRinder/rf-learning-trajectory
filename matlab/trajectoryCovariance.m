function K = trajectoryCovariance(thetai, thetaj, trajectory, hyper)
    if isempty(thetai) %compute K from prior
        K = zeros(size(thetaj,1));
        for i = 1:size(thetaj,1)
            %disp(i);
            for j = i:size(thetaj,1)
                if i==j
                    K(i,j) = 0.5; %adding the transposed matrix later sums the diagonal elements to 1
                else
                    K(i,j) = kernelTrajectory(thetaj(i,:), thetaj(j,:), trajectory(i,:), trajectory(j,:), hyper);
                end
            end
        end
        K = K + K'; %K is symmetric because of symmetric kernel
    else %compute K from prior and posterior
        K = zeros(1,size(thetaj,1));
        for j = 1:size(thetaj,1)
            K(1,j) = kernelTrajectory(thetai, thetaj(j,:), [], trajectory(j,:), hyper);
        end
    end
end