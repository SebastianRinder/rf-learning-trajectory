
for i=1:size(ret,1)
    O(:,i) = ret{i,1}.objective;
end

plot(mean(O,2));
