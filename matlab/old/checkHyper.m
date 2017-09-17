function checkHyper()
    load('ret2.75_3.25.mat','ret');
    
    for i = 1:size(ret,1)
        if ~isempty(ret{i,1})
            x(i,1) = ret{i,1}.hyper;
            y(i,1) = sum(ret{i,1}.eta);
        else
            break;
        end
    end
    plot(x,y,'.');
end

