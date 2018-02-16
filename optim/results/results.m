%close all;
figure
hold on;
addpath('../shadedErrorBar');
addpath('../figureSize');

kernel = {'sexp','matern52','trajectory'};
lineColor = {'-g', '-b', '-r'};

for kernelIdx = 1:3
    Y = [];
    for i=1:32
        loadStr = sprintf('%s_%02d.mat', kernel{kernelIdx} ,i);
        try
            load(loadStr);
            Y(:,i) = ret.knownY;
        catch me
            break;
        end        
    end
    if i > 1
        ymean = movmean(mean(Y,2),5);
        ystd = movmean(std(Y,0,2),5)/5;
        l(kernelIdx) = shadedErrorBar(1:size(Y,1),ymean,ystd,'lineprops',{lineColor{kernelIdx}, 'LineWidth',2},'transparent',1);
        %l(kernelIdx) = shadedErrorBar(1:size(Y,1),mean(Y,2),std(Y,0,2),'lineprops',lineColor{kernelIdx},'transparent',1);
    end
end
grid on;
axis([0, size(ymean,1), ylim])
legend([l(3).mainLine l(2).mainLine l(1).mainLine],'trajectory kernel','matern 5/2 kernel','squared exponential kernel');

xlabel('number of evaluations')
ylabel('cumulative reward')

figuresize(13,10,'cm')
