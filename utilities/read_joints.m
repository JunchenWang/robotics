function out = read_joints(port)


data = read(port, 1, "uint8");
while port.NumDatagramsAvailable > 0
     data = read(port, 1, "uint8"); 
end
%% ===== bytes -> JSON string =====
json_str = native2unicode(data.Data, 'UTF-8');
out = unpack_joint_state(json_str);
end
