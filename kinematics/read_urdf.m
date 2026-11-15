function robot = read_urdf(filename,flange_link, base_link)

if nargin < 2
    robot = get_robot_chain(importrobot(filename));
elseif nargin < 3
    robot = get_robot_chain(importrobot(filename), flange_link);
else
    robot = get_robot_chain(importrobot(filename), flange_link, base_link);
end

end