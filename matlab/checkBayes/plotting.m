function plotting(observedX, acqFcn, gprMdl, minFMean, myFig)
    figure(myFig);
    set(myFig, 'Name', num2str(size(observedX,1)));
    clf;
    
    load('gaussMix.mat');
    
    EI = zeros(100,100);
    meanY = zeros(100,100);
    var = zeros(100,100);
    for j = 1:100
        for k = 1:100
            [EI(j,k), meanY(j,k), var(j,k)] = acqFcn([X1(j,k), X2(j,k)], gprMdl, minFMean);
        end
    end
    EI = -EI;

    subplot(2,2,1);
    contour(X1,X2,obj);
    hold on;
    plot(observedX(:,1),observedX(:,2),'r+');
    plot(observedX(end,1),observedX(end,2),'k*');
    colorbar;

    subplot(2,2,2);
    contour(X1,X2,meanY);
    colorbar;

    subplot(2,2,3);
    contour(X1,X2,EI);
    hold on;
    plot(observedX(:,1),observedX(:,2),'r+');
    plot(observedX(end,1),observedX(end,2),'k*');
    colorbar;

    subplot(2,2,4);
    surf(X1,X2,EI);
    hold on;
    z = -acqFcn(observedX(end,:), gprMdl, minFMean);
    plot3(observedX(end,1),observedX(end,2), z, 'ro', 'LineWidth', 5);
    colorbar;
    
    keyboard;
end

