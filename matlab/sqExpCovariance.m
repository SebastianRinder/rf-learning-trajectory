function K = sqExpCovariance(Xm, Xn, opts)
    try
        K = opts.hyper(1).*exp((-0.5./opts.hyper(2)).*(pdist2(Xm,Xn).^2));
    catch me
        keyboard;
    end
end