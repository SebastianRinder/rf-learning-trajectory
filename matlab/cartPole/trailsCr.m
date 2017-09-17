% load('ret.mat','ret');
% 
% for trial = 1:size(ret,1)
%     E(:,trial) = ret{trial,1}.eta;
% end
% disp();
% 


for i=1:100
    a(i) = A{i,1}.action;
end