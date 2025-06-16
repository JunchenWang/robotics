function [T, f, h, y, R, Q, Cx, Cy, N] = FilterInputPreprocess(f, h, y, R, Q, Cx, Cy, x0)
if iscell(y)
    T = length(y);
else   
    [M T] = size(y);
    y = mat2cell(y, M, ones(1, T));
end
if nargin == 8
    N = length(x0);
end
if ~iscell(f)
    fs = cell(1, T);
    for i = 1 : T
        fs{i} = f;
    end
    f = fs;
end
if ~iscell(h)
    hs = cell(1, T);
    for i = 1 : T
        hs{i} = h;
    end
    h = hs;
end
if ~iscell(h)
    hs = cell(1, T);
    for i = 1 : T
        hs{i} = h;
    end
    h = hs;
end
if ~iscell(R)
    R = R(:, :);
    [DimR Num] = size(R);
    if Num == DimR
        R = repmat(R, 1, T);
    end
    R = mat2cell(R, DimR, DimR * ones(1, T));
end

if ~iscell(Q)
    Q = Q(:, :);
    [DimQ Num] = size(Q);
    if Num == DimQ
        Q = repmat(Q, 1, T);
    end
    Q = mat2cell(Q, DimQ, DimQ * ones(1, T));
end

if ~iscell(Cx)
    if isempty(Cx)
        Cx = cell(1, T);
    else
        Cx = Cx(:, :);
        [len Num] = size(Cx);
        if Num == 1
            Cx = repmat(Cx, 1, T);
        end
        Cx = mat2cell(Cx, len, ones(1, T));
    end
end

if ~iscell(Cy)
 if isempty(Cy)
        Cy = cell(1, T);
    else
        Cy = Cy(:, :);
        [len Num] = size(Cy);
        if Num == 1
            Cy = repmat(Cy, 1, T);
        end
        Cy = mat2cell(Cy, len, ones(1, T));
 end
end
