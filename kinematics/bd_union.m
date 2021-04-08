function bounds = bd_union(bds)
    n = length(bds)/2;
    if n == 1
        bounds = bds;
        return;
    end
    bounds = zeros(1,2*n);
    for i = 1 : n - 1
        for j = 1 : n - i
            if bds(2*j-1) > bds(2*j+1)
                a = bds(2*j-1);
                b = bds(2*j);
                bds(2*j-1) = bds(2*j+1);
                bds(2*j) = bds(2*j+2);
                bds(2*j+1) = a;
                bds(2*j+2) = b;
            end
        end
    end
    a = bds(1);
    b = bds(2);
    cnt = 0;
    for i = 2 : n
        if b < bds(2*i-1)
            cnt = cnt + 1;
            bounds(2*cnt - 1) = a;
            bounds(2*cnt) = b;
            a = bds(2*i-1);
            b = bds(2*i);
        else
            b = bds(2*i);
        end
    end
     cnt = cnt + 1;
     bounds(2*cnt - 1) = a;
     bounds(2*cnt) = b;
     bounds(2*cnt + 1:end) = [];
end