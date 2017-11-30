function K = scaleKernel(D, opts)
    K = opts.hyper(:,1) .* exp(-D./opts.hyper(:,2));
end

