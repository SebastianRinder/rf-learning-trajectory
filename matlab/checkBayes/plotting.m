function plotting(knownX, knownY, opts)
    selectFigure('check Bayes');
    clf;
    
    oldKnownX = knownX(1:end-1,:);
    oldKnownY = knownY(1:end-1,1);
    
    load('gaussMix.mat');
    
    EI = zeros(100,100);
    meanY = zeros(100,100);
    var = zeros(100,100);
    for j = 1:100
        for k = 1:100
            [EI(j,k), meanY(j,k), var(j,k)] = expectedImprovement([X1(j,k), X2(j,k)], oldKnownX, oldKnownY, opts);
        end
    end
    
    subplot(2,2,1);
    surf(X1,X2,obj);
%     contour(X1,X2,obj);
    hold on;
    plot3(oldKnownX(:,1),oldKnownX(:,2),oldKnownY(:,1),'r+');
    plot3(knownX(end,1),knownX(end,2),knownY(end,1),'g*');
    colorbar;

    subplot(2,2,2);
    surf(X1,X2,meanY);
    hold on;
    plot3(oldKnownX(:,1),oldKnownX(:,2),oldKnownY(:,1),'r+');
    plot3(knownX(end,1),knownX(end,2),knownY(end,1),'g*');
%     contour(X1,X2,meanY);
    colorbar;

    subplot(2,2,3);
    surf(X1,X2,var);
    hold on;
    plot3(oldKnownX(:,1),oldKnownX(:,2),oldKnownY(:,1),'r+');
    plot3(knownX(end,1),knownX(end,2),knownY(end,1),'g*');
    
%     contour(X1,X2,EI);
%     hold on;
%     plot(oldKnownX(:,1),oldKnownX(:,2),'r+');
%     plot(knownX(end,1),knownX(end,2),'k*');
    colorbar;

    subplot(2,2,4);
    surf(X1,X2,EI);
    hold on;
    z = expectedImprovement(knownX(end,:), oldKnownX, oldKnownY, opts);
    plot3(knownX(end,1),knownX(end,2), z, 'ro', 'LineWidth', 5);
    colorbar;
    
    selectFigure('EI vals');
    clf;
    plot(sort(EI(:))); 
    keyboard;
end

