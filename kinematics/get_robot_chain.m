function [robot, joint_names] = get_robot_chain(robot_tree, flange_link, base_link)
% 创建从base_link 到 flange_link的运动链
robot = [];
joint_names = {};
%% find base and flange
if nargin < 3 || strcmp(base_link, robot_tree.BaseName)
    base = robot_tree.Base;
else
    base_id = find(strcmp(robot_tree.BodyNames, base_link));
    if isempty(base_id)
        disp('not a valid base');
        return;
    else
        base = robot_tree.Bodies{base_id};
    end
end

if nargin < 2 
    flange = robot_tree.Bodies{end};
else
    flange_id = find(strcmp(robot_tree.BodyNames, flange_link));
    if isempty(flange_id)
        disp('not a valid flange');
        return;
    else
        flange = robot_tree.Bodies{flange_id};
    end
end

%% check chain
link = flange;
flag = 0;
maxN = length(robot_tree.Bodies);
body_list = cell(1, maxN);
head = maxN + 1;
while(~isempty(link))
    if(link == base)
        flag = 1;
        break;
    end 
    head = head - 1;
    body_list{head} = link;
    link = link.Parent;
end
body_list = body_list(head:end);
n = length(body_list);
if ~flag || n == 0
    disp('not a valid chain');
    return;
end

[robot, joint_names] = build_robot_chain(body_list);
end
