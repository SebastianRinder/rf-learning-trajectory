function hyper = setHypers(yTrain, lb, ub, isDeterministic)
    hyper.l = 1e-6*mean((ub - lb) / 2);
    hyper.f = std(yTrain)/sqrt(2);
    if hyper.f == 0
        hyper.f = 1;
    end

    if isDeterministic
        hyper.noise = max(1e-8, std(yTrain)*1e-4);
        hyper.noiseLB = max(hyper.noise, 1e-6);
    else
        %%To optimize
        hyper.noiseLB = max(1e-8, std(yTrain)*1e-2);
        hyper.noise = max(hyper.noiseLB, 0.1*std(yTrain)/5);
        hyper.noiseLB = max(hyper.noiseLB, 1e-6);
    end
end

