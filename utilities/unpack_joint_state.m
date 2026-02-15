function out = unpack_joint_state(jsonStr)

js = jsondecode(jsonStr);

out = struct();
out.flag  = int32(js.flag);
out.strip = int32(js.strip);
out.stamp.s = int64(js.stamp.s);
out.stamp.ns = int32(js.stamp.ns);
if out.strip <= 0
    error('strip 非法');
end

if isstring(js.names)
    out.names = cellstr(js.names);
else
    out.names = js.names;
end

dims = double(js.dims(:).');
pos = double(js.pos(:).');

if numel(out.names) ~= numel(dims)
    error('names 与 dims 长度不一致');
end

expected = sum(dims) * out.strip;
if expected ~= numel(pos)
    error('sum(dims)*strip 与 pos 长度不一致');
end

N = numel(dims);
out.values = cell(1, N);
k = 1;
for i = 1:N
    len = dims(i) * out.strip;
    out.values{i} = reshape(pos(k:k+len-1), dims(i), out.strip);
    k = k + len;
end
out.dims = dims;
end
