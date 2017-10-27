function K = sqExpCovariance(Xm, Xn, opts)
    hyper = opts.hyper;
    K = hyper.f.*exp((-0.5.*hyper.l).*(pdist2(Xm,Xn).^2));
end