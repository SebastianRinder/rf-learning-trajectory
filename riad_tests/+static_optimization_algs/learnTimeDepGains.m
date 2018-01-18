%load proTrials
clear variables;
infname = 'proTrialsMatSupraLowVar';
load (infname);

firstTimeStep = 1;
maxTimeStep = 520;
lengthSupraStep = 10;
dof = size(proTrials.trials{1}.jointPos, 2);
onlyUseBias = 1;
stateDim = 1 + 2 * dof + 3;

%cleaning up trials with ball_reset in the middle of the trajectory (why does it happen?)
for j = 1:length(proTrials.trials)
    if(proTrials.trials{j}.hitPointIndex > maxTimeStep)
        proTrials.trials{j} = {};
    end
end
proTrials.trials = proTrials.trials(~cellfun('isempty', proTrials.trials));

nbSamples = length(proTrials.trials);
nbSamples = 1;
regularization = .0000060;
gains = cell(ceil((maxTimeStep - firstTimeStep + 1) / lengthSupraStep), 1);
cholC = cell(ceil((maxTimeStep - firstTimeStep + 1) / lengthSupraStep), 1);
proTrials.supraSteps = cell(maxTimeStep - firstTimeStep + 1, 1);
unexplainedVariance = zeros(maxTimeStep - firstTimeStep + 1, dof);
currentSupraStep = 0;
proTrials.reuseGainsMode = 0;

    
maxAbsGains = zeros(stateDim, dof);
%------------------------------ Reuse gains for same supra step
if(proTrials.reuseGainsMode)
    for t = firstTimeStep:lengthSupraStep:maxTimeStep
        currentSupraStep = currentSupraStep + 1;
        lengthT = min(lengthSupraStep, maxTimeStep - t + 1);
        X = zeros(lengthT * nbSamples, stateDim);
        U = zeros(lengthT * nbSamples, dof);
        for tt = 0:(lengthT - 1)
            for j = 1:nbSamples
                X(j + tt * nbSamples, :) = [1 proTrials.trials{j}.jointPos(t + tt, :) ...
                    proTrials.trials{j}.jointVel(t + tt, :) ...
                    proTrials.trials{j}.ballTrajectory(t + tt, :)];
                U(j + tt * nbSamples, :) = proTrials.trials{j}.u(t + tt + 1, :);
            end
        end
        [gains{currentSupraStep}, cholC{currentSupraStep}] = static_optimization_algs.Normal.wmleLinearMean(X, U, [], regularization);
        maxAbsGains = max(maxAbsGains, gains{currentSupraStep});
        cholC{currentSupraStep} = chol(cholC{currentSupraStep});
        for tt = 0:(lengthT - 1)
            unexplainedVariance(t + tt - firstTimeStep + 1, :) = var(X((1 + tt * nbSamples):(nbSamples + tt * nbSamples), :)...
                * gains{currentSupraStep} - U((1 + tt * nbSamples):(nbSamples + tt * nbSamples), :))...
                ./ var(U((1 + tt * nbSamples):(nbSamples + tt * nbSamples), :));
            proTrials.supraSteps{t + tt - firstTimeStep + 1} = currentSupraStep;
        end
    end
else % Reuse actions for same supra steps
    for t = firstTimeStep:lengthSupraStep:maxTimeStep
        currentSupraStep = currentSupraStep + 1;
        lengthT = min(lengthSupraStep, maxTimeStep - t + 1);
        X = zeros(nbSamples, stateDim);
        U = zeros(lengthT, nbSamples, dof);
        for j = 1:nbSamples
            X(j, :) = [1 proTrials.trials{j}.jointPos(t, :) ...
                proTrials.trials{j}.jointVel(t, :) ...
                proTrials.trials{j}.ballTrajectory(t, :)];
            for tt = 1:lengthT
                U(tt, j, :) = proTrials.trials{j}.u(t + tt, :);
            end
        end
        meanU = permute(mean(U, 1), [2 3 1]);
        if(onlyUseBias)
            gains{currentSupraStep} = [meanU; zeros(stateDim - 1, dof)];
            cholC{currentSupraStep} = eye(dof);
        else
            [gains{currentSupraStep}, cholC{currentSupraStep}] = static_optimization_algs.Normal.wmleLinearMean(X, meanU, [], regularization);
        end
        maxAbsGains = max(maxAbsGains, gains{currentSupraStep});
        cholC{currentSupraStep} = chol(cholC{currentSupraStep});
        for tt = 1:lengthT
            if(tt == 1)
                prediction = X * gains{currentSupraStep};
            end
            unexplainedVariance(t + tt - firstTimeStep, :) = var(prediction - permute(U(tt, :, :), [2 3 1]))...
                ./ var(permute(U(tt, :, :), [2 3 1]));
            proTrials.supraSteps{t + tt - firstTimeStep} = currentSupraStep;
        end
    end
end

proTrials.gains = gains;
proTrials.cholC = cholC;
proTrials.unexplainedVariance = unexplainedVariance;
proTrials.firstTimeStep = firstTimeStep;
proTrials.maxTimeStep = maxTimeStep;
proTrials.trials = {};
outfname = [infname 'OneShotWithGains' num2str(lengthSupraStep)];
if(onlyUseBias)
    outfname = [outfname '_biasOnly'];
end
if(firstTimeStep == 1)
    outfname = [outfname '_noWait'];
end
save(outfname, 'proTrials');

%------------------------------ 
% X = [];
% U = [];
% supraStep = 1;
% supraSteps = cell(maxTimeStep, 1);
% for t = 1:maxTimeStep
%     t
%     for j = 1:nbSamples %add every trajectory of current time step
%         X = [X; 1 proTrials.trials{j}.jointPos(t, :) ...
%             proTrials.trials{j}.jointVel(t, :) proTrials.trials{j}.ballTrajectory(t, :)];
%         U = [U; proTrials.trials{j}.u(t+1, :)];
%     end 
%     cholM = chol(X' * X + regularization * eye(size(X,2)));
%     K = (cholM \ (X / cholM)') * U;
%     unexplainedVariance(t, :) = var(X(end-nbSamples+1:end, :) * K...
%         - U(end-nbSamples+1:end, :)) ./ var(U(end-nbSamples+1:end, :));
%     maxVar = max(unexplainedVariance(t, :))
%     if((maxVar > .3) && (size(X, 1) > nbSamples))
%         warning('unex var larger than thresh');
%         X = X(end-nbSamples+1:end, :);
%         U = U(end-nbSamples+1:end, :);
%         cholM = chol(X' * X + regularization * eye(size(X,2)));
%         K = (cholM \ (X / cholM)') * U;
%         gains{supraStep} = K;
%         unexplainedVariance(t, :) = var(X * K - U) ./ var(U);
%         supraStep = supraStep + 1;
%     end
%     supraSteps{t} = supraStep;
% end
% gains{supraStep} = K;
% gains = gains(~cellfun('isempty',gains));
% proTrials.gains = gains;
% proTrials.covs = covs;
% proTrials.unexplainedVariance = unexplainedVariance;
% proTrials.supraSteps = supraSteps;
% save proTrialsMatSupraLowVarWithGains proTrials

