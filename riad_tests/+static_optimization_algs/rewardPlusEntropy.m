function [res, reward, entropy] = rewardPlusEntropy(beta, sigma)
reward = normcdf(1, 0, sigma) - normcdf(-1, 0, sigma);
reward = 2 * reward - 1;
entropy = .5 * log(2 * pi * exp(1) * sigma^2);
res = reward + beta * entropy;
end

