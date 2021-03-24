function angles = inverse_kin_kuka(R, t, cfg)
% inverse_kin_kuka kuka med的运动学逆解,冗余信息为kesai
% cfg signs of axis 2,4,6
% kesai redundancy
lower = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
upper = -lower;
eps1 = 1e-6;
eps2 = 1e-6;
z = [0, 0, 1]';
d1 = 340;
d3 = 400;
d5 = 400;
d7 =126;
p02 = [0, 0, d1]';
p67 = [0, 0, d7]';
p26 = t - p02 - R * p67;
l_p26 = norm(p26);
p26_hat = p26 / l_p26;
l2_p26 = p26'*p26;
theta3 = 0;
theta4 = cfg(2) * acos((l2_p26 - d3*d3 - d5*d5) / (2*d3*d5));
R34 = [cos(theta4), 0, -sin(theta4);0, 1, 0;sin(theta4), 0, cos(theta4)];
if norm(cross(p26_hat, z)) < eps1
    theta1 = 0;
else
    theta1 = atan2(p26(2), p26(1));
end
phi = real(acos((d3*d3 + l2_p26 - d5*d5)/(2*d3*l_p26)));
theta2 = atan2(sqrt(p26(1)^2 + p26(2)^2), p26(3)) + cfg(2)*phi;
T03 = forward_kin_kuka([theta1, theta2, theta3]);
R03 = T03(1:3, 1:3);
% shoulder
As = so_w(p26_hat) * R03;
Bs = -so_w(p26_hat) * As;
Cs = p26_hat * p26_hat' * R03;
% wrist
Aw = R34' * As' * R;
Bw = R34' * Bs' * R;
Cw = R34' * Cs' * R;

bounds = anlysis_kesai(As, Bs, Cs, Aw, Bw, Cw, cfg, lower, upper);
if isempty(bounds)
    error('can not find kesai!');
end
kesai = (bounds(1) + bounds(2)) / 2;
theta2 = cfg(1) * real(acos(As(3,3) * sin(kesai) + Bs(3,3) * cos(kesai) + Cs(3,3)));
if abs(theta2) < eps2
    % theta1 + theta3;
    theta1 = atan2(As(2,1) * sin(kesai) + Bs(2,1) * cos(kesai) + Cs(2,1),...
        As(1,1) * sin(kesai) + Bs(1,1) * cos(kesai) + Cs(1,1));
    theta3 = 0;
else
    theta1 = atan2(cfg(1) * (As(2,3) * sin(kesai) + Bs(2,3) * cos(kesai) + Cs(2,3)),...
        cfg(1) * (As(1,3) * sin(kesai) + Bs(1,3) * cos(kesai) + Cs(1,3)));
    theta3 = atan2(cfg(1) * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
        -cfg(1) * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
end
theta6 = cfg(3) * real(acos(Aw(3,3) * sin(kesai) + Bw(3,3) * cos(kesai) + Cw(3,3)));
if abs(theta6) < eps2
    % theta5 + theta7 ;
    theta7 = 0;
    theta5 = atan2(Aw(2,1) * sin(kesai) + Bw(2,1) * cos(kesai) + Cw(2,1),...
        Aw(1,1) * sin(kesai) + Bw(1,1) * cos(kesai) + Cw(1,1));
else
    theta5 = atan2(cfg(3) * (Aw(2,3) * sin(kesai) + Bw(2,3) * cos(kesai) + Cw(2,3)),...
        cfg(3) * (Aw(1,3) * sin(kesai) + Bw(1,3) * cos(kesai) + Cw(1,3)));
    theta7 = atan2(cfg(3) * (Aw(3,2) * sin(kesai) + Bw(3,2) * cos(kesai) + Cw(3,2)),...
        -cfg(3) * (Aw(3,1) * sin(kesai) + Bw(3,1) * cos(kesai) + Cw(3,1)));
end
angles = [theta1, theta2, theta3, theta4, theta5, theta6, theta7] / pi * 180;
end


function theta = kesai2theta_hj(kesai, a, b, c, cfg)
    theta = cfg*real(acos(a*sin(kesai) + b*cos(kesai) + c));
end

function kesai = theta2kesai_hj(theta, a, b, c)
    delta = a^2+b^2-(c-cos(theta))^2;
    kesai = [2*atan((a-sqrt(delta))/(cos(theta)+b-c)), 2*atan((a+sqrt(delta))/(cos(theta)+b-c))];
    if kesai(1) > kesai(2)
        tem = kesai(1);
        kesai(1) = kesai(2);
        kesai(2) = tem;
    end
end

function bounds = kesai_range_hj(str, a, b, c, cfg, lower, upper)
    y = @(kesai) kesai2theta_hj(kesai, a, b, c, cfg);
    k = @(theta) theta2kesai_hj(theta, a, b, c);
    x = linspace(-pi, pi, 100);
    x = x(2:end);
    figure;
    plot(x, y(x));
    title(str);
    grid on;
    delta = a^2+b^2;
    if abs(delta) < 1e-12 % a=0,b=0 theta is appx. constant wrt. kesasi
        theta = y(0);
        if theta > lower && theta < upper
            bounds =[-pi, pi];
        else
            bounds = [];
        end
    else % two stationary points
        kesai = atan(a/b);
        theta = y(kesai);
        theta2 = y(kesai + pi);
        if abs(theta) < 1e-6 % singular position
            bounds=[-pi, theta - 0.1, theta + 0.1, pi];
            if cfg > 0 && theta2 > upper
                ks = k(upper);
                bounds(1) = ks(1);
                bounds(4) = ks(2);
            elseif cfg < 0 && theta2 < lower
                ks = k(lower);
                bounds(1) = ks(1);
                bounds(4) = ks(2);
            end
        else
            d2 = (b*cos(kesai)+a*sin(kesai))/sin(theta);
            if d2 < 0 %maxima
               if cfg > 0 && theta2 >= upper || cfg < 0 && theta <= lower
                   bounds = []; 
               elseif cfg > 0 && theta > upper
                   ks = k(upper);
                   bounds(1) = -pi;
                   bounds(2) = ks(1);
                   bounds(3) = ks(2);
                   bounds(4) = pi;
               elseif cfg < 0 && theta2 < lower
                   bounds = k(lower);
               else
                   bounds = [-pi, pi];
               end 
            else % minima
               if cfg > 0 && theta >= upper || cfg < 0 && theta2 <= lower
                   bounds = [];  
               elseif cfg > 0 && theta2 > upper
                   bounds = k(upper);
               elseif cfg < 0 && theta < lower
                   ks = k(lower);
                   bounds(1) = -pi;
                   bounds(2) = ks(1);
                   bounds(3) = ks(2);
                   bounds(4) = pi;
               else
                   bounds = [-pi, pi];
               end 
            end    
        end
    end
end


function theta = kesai2theta_pj(kesai, an, bn, cn, ad, bd, cd, cfg)
    theta = atan2(cfg*(an*sin(kesai)+bn*cos(kesai)+cn),...
                       cfg*(ad*sin(kesai)+bd*cos(kesai)+cd));
end

function kesai = theta2kesai_pj(theta, an, bn, cn, ad, bd, cd)
    ap = (cd-bd)*tan(theta)+(bn-cn);
    bp = 2*(ad*tan(theta)-an);
    cp = (bd+cd)*tan(theta)-(bn+cn);
    delta2 = bp^2 - 4*ap*cp;
    kesai = [2*atan((-bp-sqrt(delta2)) / (2*ap)), 2*atan((-bp+sqrt(delta2)) / (2*ap))];
    if kesai(1) > kesai(2)
        tem = kesai(1);
        kesai(1) = kesai(2);
        kesai(2) = tem;
    end
end


function bounds = bd_intersection(bd1, bd2)
    if isempty(bd1) || isempty(bd2)
        bounds = [];
        return;
    end
    bounds = zeros(1,4);
    cnt = 1;
    for i = 1 : 2 : length(bd1)
        for j = 1 : 2 : length(bd2)
            a = max(bd1(i), bd2(j));
            b = min(bd1(i+1), bd2(j+1));
            if a<=b
                bounds(cnt) = a;
                bounds(cnt+1) = b;
                cnt = cnt + 2;
            end
        end
    end
    bounds(cnt:end) = [];
end

function bounds = kesai_range_pj(str, an, bn, cn, ad, bd, cd, cfg, lower, upper)
% bn = -bn;
% bd = -bd;
% cn = bn +1e-6;
% cd = bd;
    y = @(kesai) kesai2theta_pj(kesai, an, bn, cn, ad, bd, cd, cfg);
    k = @(theta) theta2kesai_pj(theta, an, bn, cn, ad, bd, cd);
    x = linspace(-pi, pi, 1000);
    x = x(2:end);
    figure;
    plot(x, y(x));
    title(str);
    grid on;
    % theta = atan2(u, v)
    at = cn*bd-bn*cd; %原论文有错，不带cfg
    bt = an*cd-cn*ad;
    ct = an*bd-bn*ad;
    if abs(at) < 1e-6 && abs(bt) < 1e-6 && abs(ct) < 1e-6 % theta is constant wrt. kesai
        theta = y(0);
        if upper >= theta && lower <= theta
            bounds = [-pi, pi];
            return;
        else
            bounds=[];
            return;
        end
    end
    delta = at^2 + bt^2 - ct^2;
    if abs(delta) <= 1e-6
        tol = 0.1; % 5.7 degree
        if bt ~= ct
            kesai_s = 2*atan(at/(bt-ct)); % singular arm angle
        else
            error('check it 1!');
%             kesai_s = 2*atan(-(bt+ct)/(2*at));
        end
        kesai1 = kesai_s - tol;
        kesai2 = kesai_s + tol;
        bounds(1) = -pi;
        bounds(2) = kesai1;
        bounds(3) = kesai2;
        bounds(4) = pi;
        theta1 = y(kesai1);
        d1 = at*sin(kesai1)+bt*cos(kesai1)+ct;
        theta2 = y(kesai2);
%         d2 = at*sin(kesai2)+bt*cos(kesai2)+ct;  
        if theta1 > theta2 % jump -pi
            if d1 > 0 % 1
                if theta1 > upper
                    for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-3
                            break;
                        end
                    end
                    bounds(2) = kesai;
                end
                if theta2 < lower
                    for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-3
                            break;
                        end
                    end
                    bounds(3) = kesai;
                end
            else % 2
                flag1 = 0; flag2 = 0;
                if theta1 < upper
                    for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-3
                            break;
                        end
                    end
                    bounds(1) = kesai;
                else
                    flag1 = 1;
                end
                if theta2 > lower
                    for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-3
                            break;
                        end
                    end
                    bounds(4) = kesai;
                else
                    flag2=1;
                end
                if flag2 == 1
                    bounds(3:4) = [];
                end
                if flag1 == 1
                    bounds(1:2) = [];
                end
            end
        else % jump +pi
           if d1 < 0 % 1
                if theta1 < lower
                    for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-3
                            break;
                        end
                    end
                    bounds(2) = kesai;
                end
                if theta2 > upper
                    for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-3
                            break;
                        end
                    end
                    bounds(3) = kesai;
                end
            else % 2
                 flag1 = 0; flag2 = 0;
                if theta1 > lower
                    for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-3
                            break;
                        end
                    end
                    if kesai < kesai1
                        bounds(1) = kesai;
                    end
                else
                    flag1 = 1;
                end
                if theta2 < upper
                    for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-3
                            break;
                        end
                    end
                    if kesai > kesai2
                        bounds(4) = kesai;
                    end
                else
                    flag2 = 1;
                end
                 if flag2 == 1
                    bounds(3:4) = [];
                end
                if flag1 == 1
                    bounds(1:2) = [];
                end
            end
        end  
    elseif delta > 0
%           kesai1 = asin(-ct/(sqrt(at^2+bt^2))) - atan2(bt, at);
%           kesai2 = pi - asin(-ct/(sqrt(at^2+bt^2))) - atan2(bt, at);
        if bt ~= ct
            kesai1 = 2*atan((at + sqrt(delta)) / (bt - ct));
            kesai2 = 2*atan((at - sqrt(delta)) / (bt - ct));
        else
            error('bt == ct!!!');
        end
        if at*cos(kesai1)-bt*sin(kesai1) < 0 % local maxmia
            tem = kesai2;
            kesai2 = kesai1;
            kesai1 = tem;
        end
        % now kesai1: local minima, kesai2: local maxmia
        y1 = y(kesai1);
        y2 = y(kesai2);
        if y1 < y2 % no jump
            bounds = [-pi, pi];
            if y1 < lower
                kesai = k(lower);
                bounds(1) = -pi;
                bounds(2) = kesai(1);
                bounds(3) = kesai(2);
                bounds(4) = pi;
            end
            if y2 > upper
                kesai = k(upper);
                if length(bounds) == 4
                    bounds(2) = min(bounds(2),kesai(1));
                    bounds(3) = max(bounds(3), kesai(2));
                else
                    bounds(1) = -pi;
                    bounds(2) = kesai(1);
                    bounds(3) = kesai(2);
                    bounds(4) = pi;
                end
            end
        else % with jump
            bounds =[];
            if y1 < upper
                kesai = k(upper);
                bounds(1) = kesai(1);
                bounds(2) = kesai(2);
            end 
            if y2 > lower
                kesai = k(lower);
                bounds(3) = kesai(1);
                bounds(4) = kesai(2);  
            end
        end
    else % delta < 0
        for kesai = k(lower)
             if abs(y(kesai) - lower) < 1e-3
                 break;
             end
        end
        kesai1 = kesai;
       for kesai = k(upper)
             if abs(y(kesai) - upper) < 1e-3
                 break;
             end
       end
        kesai2 = kesai;
        if kesai1 > kesai2
            tem = kesai1;
            kesai1 = kesai2;
            kesai2 = tem;
        end
        
        theta = y((kesai1 + kesai2) / 2);
        if theta > lower && theta < upper % no jump
            bounds(1) = kesai1;
            bounds(2) = kesai2;
        else % with jump
            bounds(1) = -pi;
            bounds(2) = kesai1;
            bounds(3) = kesai2;
            bounds(4) = pi;
        end
    end 
end

function kesais = anlysis_kesai(As, Bs, Cs, Aw, Bw, Cw, cfg, lower, upper)
% bd = bd_intersection([1,3,5,8],[2,6]);
% an = As(3,2); bn = Bs(3,2); cn = Cs(3,2);
% ad = -As(3,1); bd = -Bs(3,1); cd = -Cs(3,1);
kesais1 = kesai_range_pj('joint 1', As(2,3), Bs(2,3), Cs(2,3), As(1,3), Bs(1,3), Cs(1,3), cfg(1),lower(1), upper(1));% A1
kesais3 = kesai_range_pj( 'joint 3', As(3,2), Bs(3,2), Cs(3,2), -As(3,1), -Bs(3,1), -Cs(3,1), cfg(1),lower(3), upper(3));% A3
kesais5 = kesai_range_pj('joint 5',  Aw(2,3), Bw(2,3), Cw(2,3), Aw(1,3), Bw(1,3), Cw(1,3), cfg(2),lower(5), upper(5));% A1
kesais7 = kesai_range_pj('joint 7',  Aw(3,2), Bw(3,2), Cw(3,2), -Aw(3,1), -Bw(3,1), -Cw(3,1), cfg(3),lower(7), upper(7));% A3
kesais2 = kesai_range_hj('joint 2',  As(3,3), Bs(3,3), Cs(3,3), cfg(1),lower(2), upper(2));% A2
kesais6 = kesai_range_hj('joint 6',  Aw(3,3), Bw(3,3), Cw(3,3), cfg(3),lower(6), upper(6));% A2
kesais = bd_intersection(kesais1,kesais2);
kesais = bd_intersection(kesais,kesais3);
kesais = bd_intersection(kesais,kesais5);
kesais = bd_intersection(kesais,kesais6);
kesais = bd_intersection(kesais,kesais7);
% kesai_range_hj( As(3,3), Bs(3,3), Cs(3,3), cfg(1),lower(2), upper(2));% A2
%              u = cfg*(an*sin(kesai) + bn*cos(kesai) + cn);
%              v = cfg*(ad*sin(kesai) + bd*cos(kesai) + cd);
%             du = cfg*(an*cos(kesai) - bn*sin(kesai));
%             dv = cfg*(ad*cos(kesai) - bd*sin(kesai));
%             fenmu = u^2 + v^2;
%    
%             d = (at*sin(kesai) + bt*cos(kesai) + ct) / fenmu;
%             d2 = (at*cos(kesai)-bt*sin(kesai) - d*(2*u*du + 2*v*dv))/fenmu;
end
