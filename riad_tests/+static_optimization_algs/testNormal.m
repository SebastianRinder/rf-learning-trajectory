% rosen = static_optimization_algs.Rosenbrock;
% line = (-1:0.1:1)';
% x1 = repmat(line, length(line), 1);
% x = [x1 reshape((repmat(line, 1, length(line)))',size(x1))];
% y = rosen.f(x);
% surf(line, line, reshape(y, length(line), length(line))');
clear variables;
close all;

figure('Position', [800, 400, 800, 600]);

%covC = [1 .5; .5 1];
covC = 5 * (ones(2,2) - 2 * rand(2,2));
covC = covC' * covC;
mu = 20 * (ones(1,2) - 2 * rand(1,2));

covC2 = 5 * (ones(2,2) - 2 * rand(2,2));
covC2 = covC2' * covC2;
mu2 = 5 * (ones(1,2) - 2 * rand(1,2));

normal = static_optimization_algs.Normal(mu, covC);
normal2 = static_optimization_algs.Normal(mu2, covC2);

nbSamples = 100;
samples = normal.getSamples(nbSamples);
samples = [samples; normal2.getSamples(nbSamples)];

% Test of getLogProbas
% precision = inv(covC);
% detP = det(precision);
% normal.getLogProbas(samples) - normal.getLogProbasPrecision(samples, mu, precision, detP)

subplot(1,3,1);
hold on;
plot(samples(1:nbSamples, 1), samples(1:nbSamples, 2), '.b');
plot(samples(nbSamples+1:end, 1), samples(nbSamples+1:end, 2), '.r');
normal.plot();
normal2.plot();

subplot(1,3,2);
hold on;
plot(samples(:, 1), samples(:, 2), '.');
weightRatio = 1;
weights = [ones(nbSamples, 1); ones(nbSamples, 1) * weightRatio];
normal3 = normal.wmle(samples, weights);
normal3.plot();

subplot(1,3,3);
hold on;
plot(samples(:, 1), samples(:, 2), '.');
weightRatio = weightRatio * 10;
weights = [ones(nbSamples, 1); ones(nbSamples, 1) * weightRatio];
normal4 = normal.wmle(samples, weights);
normal4.plot();

%covC - normal.getCovariance()
