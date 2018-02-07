function D = euclideanDistance(Xm, Xn)
    m = size(Xm,1);
    n = size(Xn,1);

    D = zeros(m,n);

    if m==n && all(all(Xm == Xn)) %symmetric
        for i=1:size(Xm,1)
            D(i,i) = 0;
            k = i+1;
            for j=k:size(Xn,1)
                r = Xm(i,:) - Xn(j,:);
                D(i,j) = r*r';
            end
        end
        D = sqrt(D + D');
    else
        for i=1:size(Xm,1)
            for j=1:size(Xn,1)
                r = Xm(i,:) - Xn(j,:);
                D(i,j) = r*r';
            end
        end
        D = sqrt(D);
    end

    %pdist2(Xm,Xn); %licence error on cluster sometimes
end