%close all;
figure
hold on;
addpath('shadedErrorBar');

tmp1 = load('local_thompson_full_300_cartPole_pygym_sexp_1e-06_25-01-2018_23-29.mat')
tmp2 = load('local_thompson_full_300_cartPole_pygym_matern52_1e-06_26-01-2018_00-25.mat')

ret = tmp1.ret;
Y1 = [];
for i=1:size(ret,1)
    Y1(:,i) = ret{i,1}.knownY;
    %plot(-9:200,Y(:,i),'.b');
end
l1 = shadedErrorBar(1:size(Y1,1),mean(Y1,2),std(Y1,0,2),'lineprops','-g','transparent',1);

ret = tmp2.ret;
Y2 = [];
for i=1:size(ret,1)
    Y2(:,i) = ret{i,1}.knownY;
    %plot(-9:200,Y(:,i),'.r');
end
l2 = shadedErrorBar(1:size(Y2,1),mean(Y2,2),std(Y2,0,2),'lineprops','-b','transparent',1);

for i=1:7
    loadStr = sprintf('trajectory_%02d.mat', i);
    load(loadStr);
    Y3(:,i) = ret.knownY;
end

% ret = tmp3.ret;
% Y3 = [];
% for i=1:size(ret,1)
%     Y3(:,i) = ret{i,1}.knownY;
%     %plot(-9:200,Y3(:,i),'.g');
% end
l3 = shadedErrorBar(1:size(Y3,1),mean(Y3,2),std(Y3,0,2),'lineprops','-r','transparent',1);

legend([l3.mainLine l2.mainLine l1.mainLine],'trajectory kernel','matern 5/2 kernel','squared exponential kernel');

xlabel('bayesian optimization iteration step')
ylabel('cumulative reward')

figuresize(13,10,'cm')