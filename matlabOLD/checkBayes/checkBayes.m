function checkBayes()
    hyper = 1;
    
    policyUb = 10.*ones(1,1);
    policyLb = -1.*policyUb;

	theta(1,:) = randTheta(policyLb, policyUb);
    eta(1,1) = objective(theta(1,:));
    disp(['step : ',num2str(1) ,'  |  cumulative reward: ', num2str(eta(1,:))]);    

    r = -10:0.1:10;
    
    for i=2:100
        clf;
        K = sqExpCovariance(theta(1:i-1,:), theta(1:i-1,:), hyper);        
        toMinimize = @(x) -acquisition(x, eta(1:i-1,:), theta(1:i-1,:), @sqExpCovariance, K, hyper, []);
        theta(i,:) = globalMin(toMinimize, policyLb, policyUb);
        eta(i,:) = objective(theta(i,:));

        disp(['step : ',num2str(i) ,'  |  cumulative reward: ', num2str(eta(i,:))]);        

        k = 0;
        for j = r
            k = k + 1;
            [acq(1,k), meanY(1,k), var(1,k)] = acquisition(j, eta(1:i-1,:), theta(1:i-1,:), @sqExpCovariance, K, hyper, []);
        end
        
        %acq = -acq;
        
        subplot(2,1,1);
        errorbar(r,meanY,var,'g');
        hold on;
        plot(r,objective(r),'k');        
        plot(theta,eta,'r+');
        plot(theta(end,1),eta(end,1),'k*');
                
        subplot(2,1,2);
        plot(r,acq,'k');        
    end

end

function ret = objective(x)
    ret = sin(x);
end

