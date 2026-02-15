function json_str = pack_joint_state(flag, names, pos_cell)
% names  : cell array of joint names
% pos_cell: cell array, each element is numeric array, rows are dims

if numel(names) ~= numel(pos_cell)
    error('names 与 posCell 数量必须一致');
end

% ---- stamp（快&够用）----
t = posixtime(datetime('now','TimeZone','UTC'));
sec  = int64(floor(t));
nsec = int32(round((t - double(sec))*1e9));

% ---- dims + flat pos ----
N = numel(pos_cell);
dims = zeros(1, N, 'int32');

% 预计算总长度，避免反复扩容（关键提速点）
total = 0;
strip = size(pos_cell{1},2);
for i = 1:N
    v = double(pos_cell{i});
    sz = size(v);
    dims(i) = int32(sz(1));
    if sz(2) ~= strip
        error('strip wrong');
    end
    total = total + numel(v);
end

pos = zeros(1, total);  % flat
k = 1;
for i = 1:N
    v = double(pos_cell{i});
    m = numel(v);
    pos(k:k+m-1) = v(:).';
    k = k + m;
end
if isscalar(dims)
    dims = {dims};
end
if isscalar(pos)
    pos = {pos};
end
% ---- 组包 ----
js = struct();
js.flag = int32(flag);
js.strip = int32(strip);
js.stamp = struct('s', sec, 'ns', nsec);
js.names = names;
js.dims  = dims;
js.pos   = pos;

json_str = jsonencode(js);
end
