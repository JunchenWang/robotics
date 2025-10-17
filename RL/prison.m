
rng("default");

% 合作（抵赖）: 1, 背叛（坦白）: 2
episode_num = 200000;
% reward for player 1 
R = [-1 -10; 0 -8];
% e-greedy
eps = 0.01;
% learning rate
alpha = 0.02;
% Q值表
Q1 = zeros(2,1);
Q2 = zeros(2, 1);

total_reward1 = zeros(episode_num, 1);
total_reward2 = zeros(episode_num, 1);
action_history1 = zeros(episode_num, 1);
action_history2 = zeros(episode_num, 1);
avg_w = 20000;
avg_reward1 = zeros(episode_num, 1);
avg_reward2 = zeros(episode_num, 1);
for step = 1 : episode_num
    
    %player 1
    if rand < eps
        a1 = randi(2);
    else
        candi = find(Q1 == max(Q1));
        a1 = candi(randi(length(candi)));
    end

    %player 2
    if rand < eps
        a2 = randi(2);
    else
        candi = find(Q2 == max(Q2));
        a2 = candi(randi(length(candi)));
    end
    r1 = R(a1, a2);
    r2 = R(a2, a1);
    total_reward1(step) = r1;
    total_reward2(step) = r2;
    if(step >= avg_w)
        avg_reward1(step) = sum(total_reward1(step-avg_w + 1 : step)) / avg_w;
        avg_reward2(step) = sum(total_reward2(step-avg_w + 1 : step)) / avg_w;
    else
        avg_reward1(step) = mean(total_reward1(1:step));
        avg_reward2(step) = mean(total_reward2(1:step));
    end

    action_history1(step) = a1;
    action_history2(step) = a2;
    Q1(a1) = Q1(a1) + alpha * (r1 - Q1(a1));
    Q2(a2) = Q2(a2) + alpha * (r2 - Q2(a2));
end

p = length(find(action_history1 == 1)) / episode_num;
q = length(find(action_history2 == 1)) / episode_num;

disp('两个玩家合作的概率');
disp([p, q]);
%两个玩家每轮的平均收益
plot([avg_reward1, avg_reward2], 'LineWidth',2);
legend('player1 avg reward', 'player2 avg reward');