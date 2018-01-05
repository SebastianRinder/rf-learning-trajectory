function K = scaleKernel(D, opts)
    K = opts.hyper(:,1) .* exp(-D./opts.hyper(:,2));
    
    if opts.debugPlotting
        if size(K,1) == 10000
            selectFigure('posterior Kernel without 0 (sorted)');
            toPlot = K(K ~= 0);
            plot(sort(toPlot(:)));
            title([num2str(size(K,1)*size(K,2)), ' values in total']);
            pause(0.1);
        end
    end
end

