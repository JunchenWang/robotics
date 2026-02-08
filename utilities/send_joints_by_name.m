function send_joints_by_name(port, jointNames, jointPos, ip, no)
% UDP 发送关节角（rad），时间戳 = 秒 + 纳秒
if nargin < 5
    no = 7755;
end

if nargin < 4
    ip = '127.0.0.1';
end

 assert(numel(jointNames) == numel(jointPos), ...
        'jointNames and jointPos size mismatch');

    %% ===== 时间戳（UTC, sec + nsec）=====
    t = datetime("now", "TimeZone", "UTC");
    sec  = int64(posixtime(t));
    nsec = int64(mod(t.Second, 1) * 1e9);

    %% ===== joints 数组 =====
    N = numel(jointNames);
    joints = repmat(struct('name',"", 'pos',[]), 1, N);
    msg = struct();
    msg.stamp  = struct('sec', sec, 'nsec', nsec);
    for i = 1:N
        joints(i).name = string(jointNames{i});
        p = jointPos{i};
        % validateattributes(p, {'numeric'}, {'nonempty'}, mfilename, 'jointPos');
        % 强制转成行向量 → JSON array
        joints(i).pos = p(:).';
    end
    msg.joints = joints;

    %% ===== 发送 =====
    jsonStr = jsonencode(msg);
    write(port, uint8(jsonStr), "uint8", ip, no);
end
