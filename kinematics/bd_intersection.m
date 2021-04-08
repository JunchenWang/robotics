function bounds = bd_intersection(bd1, bd2)
if isempty(bd1) || isempty(bd2)
    bounds = [];
    return;
end
n1 = length(bd1)/2;
n2 = length(bd2)/2;
bounds = zeros(1,n1*n2);
cnt = 0;
for i = 1 : n1
    for j = 1 : n2
        a = max(bd1(2*i-1), bd2(2*j-1));
        b = min(bd1(2*i), bd2(2*j));
        if a<=b
            cnt = cnt + 1;
            bounds(2*cnt-1) = a;
            bounds(2*cnt) = b;
        end
    end
end
bounds(2*cnt+1:end) = [];
end