%close all;
figure
hold on;
addpath('../shadedErrorBar');

load('cartPole_matlab_sexp_1e-06_16-01-2018_14-46.mat')

Y1 = [];
for i=1:size(ret,1)
    Y1(:,i) = ret{i,1}.knownY;
    %plot(-9:200,Y(:,i),'.b');
end
l1 = shadedErrorBar(-9:200,mean(Y1,2),std(Y1,0,2),'lineprops','-g','transparent',1);

load('cartPole_matlab_matern52_1e-06_16-01-2018_15-08.mat')

Y2 = [];
for i=1:size(ret,1)
    Y2(:,i) = ret{i,1}.knownY;
    %plot(-9:200,Y(:,i),'.r');
end
l2 = shadedErrorBar(-9:200,mean(Y2,2),std(Y2,0,2),'lineprops','-b','transparent',1);


load('cartPole_matlab_trajectory_1e-06_16-01-2018_19-27.mat')
Y3 = [];
for i=1:size(ret,1)
    Y3(:,i) = ret{i,1}.knownY;
    %plot(-9:200,Y3(:,i),'.g');
end
l3 = shadedErrorBar(-9:200,mean(Y3,2),std(Y3,0,2),'lineprops','-r','transparent',1);

legend([l3.mainLine l2.mainLine l1.mainLine],'trajectory kernel','matern 5/2 kernel','squared exponential kernel');