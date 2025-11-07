function result = slerp(q1, q2, t)
    % 四元数球面线性插值
    % q1, q2: 输入单位四元数
    % t: 插值参数 [0, 1]
    
    % 计算点积
    dot_product = dot(q1, q2);
    
    % 确保选择最短路径
    if dot_product < 0
        q2 = -q2;
        dot_product = -dot_product;
    end
    
    
    % 计算角度
    theta = real(acos(dot_product));
    
    if abs(theta) < 1e-9
        % 如果角度很小，使用线性插值
        result = (1-t)*q1 + t*q2;
    else
        % SLERP公式
        result = (sin((1-t)*theta)/sin(theta))*q1 + (sin(t*theta)/sin(theta))*q2;
    end
    
    % 归一化
    result = result / norm(result);
end