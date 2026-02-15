function send_joints(port, joint_names, joint_pos, ip, no)

if nargin < 5
    no = 7755;
end

if nargin < 4
    ip = '127.0.0.1';
end

json_str = pack_joint_state(1, joint_names, joint_pos);
write(port, json_str, ip, no);
end
