function joints = query_joints(port, ip, no)
% joints is a row vector
if nargin < 3
    no = 7755;
end
if nargin < 2
    ip = "127.0.0.1";
end
if nargin < 1
    port = udpport("datagram");
end
json_str = pack_joint_state(-1, {'default'},{0});
write(port, uint8(json_str), "uint8", ip, no);
out = read_joints(port);
joints = out.values{1};
joints = joints(:)';
if nargin == 0
    clear port;
end
end