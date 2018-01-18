function noiseVariance = getNoiseVariance(func,lb,ub)
    func.opts.trajectoriesPerSample = 30;
    samples = randBound(lb,ub,10000);
    retVals = func.eval(samples);    
end