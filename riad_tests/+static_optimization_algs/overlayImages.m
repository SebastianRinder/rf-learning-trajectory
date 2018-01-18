function img = overlayImages(basename, minId, maxId, ext)
currI = rgb2gray(imread([basename, num2str(minId), ext])); %just for getting the size of the images
allI = zeros([(maxId-minId+1) size(currI)]);
for k = minId:maxId
    allI(k-minId+1, :, :) = rgb2gray(imread([basename, num2str(k), ext]));    
end
meanI = squeeze(mean(allI));
distI = zeros(size(allI));
for k = 1:size(allI,1)
    distI(k, :, :) = (squeeze(allI(k, :, :)) - meanI).^2;
end
[~, indM] = max(distI);
indM = squeeze(indM);
for 
    
end
end
