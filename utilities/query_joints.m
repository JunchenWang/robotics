function joints = query_joints(port)
% joints is a row vector
if nargin == 0
    port = udpport('byte');
end
writeline(port,"robot;","127.0.0.1",7755);
s = split(readline(port), ';');
joints = double(s(1:end-1))';
if nargin == 0
    clear port;
end
end