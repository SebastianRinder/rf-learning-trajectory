close all;
addpath('shadedErrorBar');

Y = [];

load('retCartpoleSexp');

for i=1:size(ret,1)
    Y(:,i) = ret{i,1}.knownY;
end

%plot(-9:40,mean(O,2));
% errorbar(-9:40,mean(O,2),std(O,0,2));
l1 = shadedErrorBar(-9:200,mean(Y,2),std(Y,0,2),'lineprops','-b','transparent',1);

% hold on;
% 
% Y = [];
% 
% load('retCartpoleTraj40_1-40');
% 
% for i=1:size(ret,1)
%     Y(:,i) = ret{i,1}.objective;
% end
% 
% % plot(-9:40,mean(O,2));
% % errorbar(-9:40,mean(O,2),std(O,0,2));
% l2 = shadedErrorBar(-9:40,mean(Y,2),std(Y,0,2),'lineprops','-r','transparent',1);

% legend([l1.mainLine l2.mainLine],'squared exponential kernel','trajectory kernel');