function K = scaleKernel(D, hyper)
    K = hyper(:,1) .* exp(-D./hyper(:,2));
end

