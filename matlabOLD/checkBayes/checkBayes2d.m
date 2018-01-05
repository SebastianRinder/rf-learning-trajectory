function checkBayes2d()
    hyper = 1;
    
    x = linspace(-2*pi,2*pi);
    y = linspace(0,4*pi);
    
    ub = [max(x), max(y)];
    lb = [min(x), min(y)];

	theta(1,:) = randTheta(lb, ub);
    eta(1,1) = objective(theta(1,:));
    disp(['step : ',num2str(1) ,'  |  cumulative reward: ', num2str(eta(1,1))]);    

    [X,Y] = meshgrid(x,y);

    close all;
    figure;
    for i=2:100
        clf;
        K = sqExpCovariance([], theta(1:i-1,:), [], hyper);        
        toMinimize = @(x) -acquisition(x, eta(1:i-1,:), theta(1:i-1,:), @sqExpCovariance, K, hyper, []);
        theta(i,:) = globalMin(toMinimize, [], lb, ub);
        eta(i,1) = objective(theta(i,:));

        disp(['step : ',num2str(i) ,'  |  cumulative reward: ', num2str(eta(i,1))]);        

        if mod(i,10) == 0
            for j = 1:100
                for k = 1:100
                    obj(j,k) = objective([X(j,k), Y(j,k)]);
                    [EI(j,k), meanY(j,k), var(j,k)] = acquisition([X(j,k), Y(j,k)], eta(1:i-1,1),...
                        theta(1:i-1,:), @sqExpCovariance, K, hyper, []);
                end
            end
%             figure;
%             histogram(K);
            %acc = -acc;

            subplot(2,2,1);
            contour(X,Y,obj);
            colorbar;

            subplot(2,2,2);
%             k = 1:2;
%             cmap = colormap;
%             m = 63/(max(eta)-min(eta));
%             for j=1:i                
%                 y = round(m * (eta(j,1) - min(eta)) + 1);
%                 plot(k,theta(j,k),'Color',cmap(y,:));
%                 hold on;
%             end
%             colorbar;
%             caxis([min(eta), max(eta)]);
            contour(X,Y,meanY);
            colorbar;

            subplot(2,2,3);
            contour(X,Y,EI);
            hold on;
            plot(theta(:,1),theta(:,2),'r+');
            plot(theta(end,1),theta(end,2),'k*');
            colorbar;
            
%             figure;
%             m = 63/(max(max(acq))-min(min(acq)));
%             for j = 1:10
%                 for k = 1:10
%                     y = round(m * (acq(j,k) - min(min(acq))) + 1);
%                     plot(1:2,[X(j,k), Y(j,k)],'Color',cmap(y,:));
%                     hold on;
%                 end
%             end
%             plot(1:2,theta(i,1:2),'k');
%             colorbar;
%             caxis([min(eta), max(eta)]);

            subplot(2,2,4);
            surf(X,Y,EI);
            colorbar;
        end
    end

end

function ret = objective(x)
    ret = sin(x(1,1))+cos(x(1,2));
    %ret = -ret;
end

