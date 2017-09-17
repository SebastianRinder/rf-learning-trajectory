load('ret.mat','ret');

for trial = 1:49 %size(ret,1)
    e = ret{trial,1}.eta;
    e = sum(e >=200);
    eta(trial,1) = e;
    
end
plot(eta);



% for i=1:100
%     a(i) = A{i,1}.action;
% end