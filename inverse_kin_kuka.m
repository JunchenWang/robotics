function [angles, bounds] = inverse_kin_kuka(R, t, cfg,lowers, uppers)
% inverse_kin_kuka kuka med的运动学逆解,自动选择kesai
% cfg signs of axis 2,4,6
eps1 = 1e-9;%not change
eps0 = 1e-6;%not change
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
theta4 = cfg(2) * real(acos((l2_p26 - d3*d3 - d5*d5) / (2*d3*d5)));
R34 = [cos(theta4), 0, -sin(theta4);0, 1, 0;sin(theta4), 0, cos(theta4)];
if abs(abs(dot(p26_hat, z)) - 1) < eps1
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

bounds = anlysis_kesai(As, Bs, Cs, Aw, Bw, Cw, cfg, lowers, uppers, eps0);
if isempty(bounds)
    angles=[];
    return;
end
kesai = choose_kesai(bounds);
theta2 = cfg(1) * real(acos(As(3,3) * sin(kesai) + Bs(3,3) * cos(kesai) + Cs(3,3)));
if abs(theta2) < eps0
    % theta1 + theta3;
    theta1and3 = atan2(As(2,1) * sin(kesai) + Bs(2,1) * cos(kesai) + Cs(2,1),...
        As(1,1) * sin(kesai) + Bs(1,1) * cos(kesai) + Cs(1,1));
    bd=bd_intersection([lowers(3), uppers(3)], theta1and3+[-uppers(1), -lowers(1)]);
    if isempty(bd)
        angles=[];
        bounds=[];
        return;
    else
        theta3=(bd(1)+bd(2))/2;
        theta1=theta1and3-theta3;
    end
else
    theta1 = atan2(cfg(1) * (As(2,3) * sin(kesai) + Bs(2,3) * cos(kesai) + Cs(2,3)),...
        cfg(1) * (As(1,3) * sin(kesai) + Bs(1,3) * cos(kesai) + Cs(1,3)));
    theta3 = atan2(cfg(1) * (As(3,2) * sin(kesai) + Bs(3,2) * cos(kesai) + Cs(3,2)),...
        -cfg(1) * (As(3,1) * sin(kesai) + Bs(3,1) * cos(kesai) + Cs(3,1)));
end
theta6 = cfg(3) * real(acos(Aw(3,3) * sin(kesai) + Bw(3,3) * cos(kesai) + Cw(3,3)));
if abs(theta6) < eps0
    % theta5 + theta7 ;
    theta5and7 = atan2(Aw(2,1) * sin(kesai) + Bw(2,1) * cos(kesai) + Cw(2,1),...
        Aw(1,1) * sin(kesai) + Bw(1,1) * cos(kesai) + Cw(1,1));
    bd=bd_intersection([lowers(7), uppers(7)], theta5and7+[-uppers(5), -lowers(5)]);
    if isempty(bd)
        angles=[];
        bounds=[];
        return;
    else
        theta7=(bd(1)+bd(2))/2;
        theta5=theta5and7-theta7;
    end
else
    theta5 = atan2(cfg(3) * (Aw(2,3) * sin(kesai) + Bw(2,3) * cos(kesai) + Cw(2,3)),...
        cfg(3) * (Aw(1,3) * sin(kesai) + Bw(1,3) * cos(kesai) + Cw(1,3)));
    theta7 = atan2(cfg(3) * (Aw(3,2) * sin(kesai) + Bw(3,2) * cos(kesai) + Cw(3,2)),...
        -cfg(3) * (Aw(3,1) * sin(kesai) + Bw(3,1) * cos(kesai) + Cw(3,1)));
end
angles = [theta1, theta2, theta3, theta4, theta5, theta6, theta7];
end

function kesai=choose_kesai(bounds)
len = bounds(2) - bounds(1);
j = 1;
for i = 3 : 2 : length(bounds)
    if bounds(i+1)-bounds(i) > len
        len = bounds(i+1)-bounds(i);
        j = i;
    end
end
kesai = (bounds(j) + bounds(j+1)) / 2;
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

function [bounds, kesai_s] = kesai_range_hj(str, a, b, c, cfg, lower, upper, eps0)
y = @(kesai) kesai2theta_hj(kesai, a, b, c, cfg);
k = @(theta) theta2kesai_hj(theta, a, b, c);
d = @(kesai) b*sin(kesai) - a*cos(kesai); % notice sign of sin(theta)
kesai_s=[];
%     x = linspace(-pi, pi, 1000);
%     figure;
%     plot(x, y(x));
%     title(str);
%     grid on;
if a==0 && b==0 % a=0,b=0 theta is constant wrt. kesasi
    theta = y(0);
    if theta >= lower && theta <= upper
        bounds =[-pi, pi];
    else
        bounds = [];
    end
else % two stationary points
    kesai = atan(a/b);
    theta = y(kesai);
    if kesai > 0
        kesai2=kesai-pi;
    else
        kesai2=kesai+pi;
    end
    theta2 = y(kesai2);
    height = abs(theta2-theta);
    if abs(theta) < abs(theta2) && abs(theta)<eps0
        kesai_s = kesai;
    elseif abs(theta2) < abs(theta) && abs(theta2)<eps0
        kesai_s = kesai2;
    end  
    if height < 1e-6
        theta = y(0);
        if theta > lower && theta < upper
            bounds =[-pi, pi];
        else
            bounds = [];
        end
        return;
    end
    [l,u] = reorder(theta, theta2);
    if upper<=l || lower>=u
        bounds=[];
        return;
    end
    bounds=[-pi, pi];
    if upper < u
        ks = k(upper);
        if d(ks(1))*sign(upper) > 0
            bounds = bd_intersection(bounds, [-pi, ks(1), ks(2), pi]);
        else
            bounds = bd_intersection(bounds, [ks(1), ks(2)]);
        end
    end
    if lower > l
        ks = k(lower);
        if d(ks(1))*sign(lower) > 0
            bounds = bd_intersection(bounds, [ks(1), ks(2)]);
        else
            bounds = bd_intersection(bounds, [-pi, ks(1), ks(2), pi]);
        end
    end
end
end


function theta = kesai2theta_pj(kesai, an, bn, cn, ad, bd, cd, cfg)
theta = atan2(cfg*(an*sin(kesai)+bn*cos(kesai)+cn),...
    cfg*(ad*sin(kesai)+bd*cos(kesai)+cd));
end

function kesai = theta2kesai_pj_2(theta, an, bn, cn, ad, bd, cd)
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

function kesai = theta2kesai_pj_1(theta, an, bn, cn, ad, bd, cd, cfg)
ap = (cd-bd)*tan(theta)+(bn-cn);
bp = 2*(ad*tan(theta)-an);
cp = (bd+cd)*tan(theta)-(bn+cn);
if abs(ap) < eps(bp)
    kesai = -cp/bp;
    return;
end
delta2 = bp^2 - 4*ap*cp;
kesai1 = 2*atan((-bp-sqrt(delta2)) / (2*ap));
kesai2 = 2*atan((-bp+sqrt(delta2)) / (2*ap));
e1 = abs(kesai2theta_pj(kesai1, an, bn, cn, ad, bd, cd, cfg) - theta);
e2 = abs(kesai2theta_pj(kesai2, an, bn, cn, ad, bd, cd, cfg) - theta);
if e1 < e2
    kesai = kesai1;
else
    kesai = kesai2;
end
end

function kesai = theta2kesai_pj_1_s(theta, an, bn, cn, ad, bd, cd, kesai_s)
ap = (cd-bd)*tan(theta)+(bn-cn);
bp = 2*(ad*tan(theta)-an);
cp = (bd+cd)*tan(theta)-(bn+cn);
delta2 = bp^2 - 4*ap*cp;
if abs(ap) < eps(bp)
    kesai = -cp/bp;
    return;
end
kesai = [2*atan((-bp-sqrt(delta2)) / (2*ap)), 2*atan((-bp+sqrt(delta2)) / (2*ap))];
if pi - abs(kesai_s) > 1e-2
    if abs(kesai(1)-kesai_s) > abs(kesai(2)-kesai_s)
        kesai = kesai(1);
    else
        kesai = kesai(2);
    end
else
    if abs(kesai(1)) < abs(kesai(2))
        kesai = kesai(1);
    else
        kesai = kesai(2);
    end
end
end

function [kesai1, kesai2] = reorder(kesai1, kesai2)
if kesai1 > kesai2
    tem = kesai1;
    kesai1 = kesai2;
    kesai2 = tem;
end
end

function bounds = bd_intersection(bd1, bd2)
if isempty(bd1) || isempty(bd2)
    bounds = [];
    return;
end
n1 = length(bd1)/2;
n2 = length(bd2)/2;
bounds = zeros(1,n1*n2);
cnt = 0;
for i = 1 : n1
    for j = 1 : n2
        a = max(bd1(2*i-1), bd2(2*j-1));
        b = min(bd1(2*i), bd2(2*j));
        if a<=b
            cnt = cnt + 1;
            bounds(2*cnt-1) = a;
            bounds(2*cnt) = b;
        end
    end
end
bounds(2*cnt+1:end) = [];
end

function bounds = kesai_range_pj(str, an, bn, cn, ad, bd, cd, cfg, lower, upper, kesai_s)
y = @(kesai) kesai2theta_pj(kesai, an, bn, cn, ad, bd, cd, cfg);
k2 = @(theta) theta2kesai_pj_2(theta, an, bn, cn, ad, bd, cd);
k1 = @(theta) theta2kesai_pj_1(theta, an, bn, cn, ad, bd, cd, cfg);
ks1 = @(theta) theta2kesai_pj_1_s(theta, an, bn, cn, ad, bd, cd, kesai_s);
tol = 0.05;
%     x = linspace(-pi, pi, 1000);
%     figure;
%     plot(x, y(x));
%     title(str);
%     grid on;
at = cn*bd-bn*cd; %原论文有错，不带cfg
bt = an*cd-cn*ad;
ct = an*bd-bn*ad;
d = @(kesai) at*sin(kesai)+bt*cos(kesai)+ct;
%     if abs(at) < eps && abs(bt) < eps && abs(ct) < eps % theta is constant wrt. kesai
if abs(at) == 0 && abs(bt)== 0 && abs(ct) == 0
    theta = y(0);
    if upper >= theta && lower <= theta
        bounds = [-pi, pi];
    else
        bounds=[];
    end
    return;
end
delta = (at^2 + bt^2 - ct^2);
delta_n = delta / (at^2 + bt^2 + ct^2);
if ~isempty(kesai_s) % singularity
    % cal tol if delta_n larger
    if delta_n>1e-3 && bt ~= ct
        kesai1 = 2*atan((at + sqrt(delta)) / (bt - ct));
        kesai2 = 2*atan((at - sqrt(delta)) / (bt - ct));
        tol = max(tol,min(abs(kesai2-kesai1), 2*pi-abs(kesai2-kesai1)));
    end
    kesai1 = kesai_s - tol;
    kesai2 = kesai_s + tol;
    theta1 = y(kesai1);
    theta2 = y(kesai2);
    [l,u] = reorder(theta1,theta2);
    if kesai1<-pi
        bounds = [kesai2, 2*pi+kesai1];
    elseif kesai2>pi
        bounds = [kesai2-2*pi, kesai1];
    else
        bounds = [-pi, kesai1, kesai2, pi];
    end
    d1 = d(kesai1);
    %         d2 = at*sin(kesai2)+bt*cos(kesai2)+ct;
    if (theta1-theta2)*d1 > 0
        if upper<=l||lower>=u
            bounds = [];
            return;
        end
        if upper>=u 
            kesai_u=kesai_s;
        else
            kesai_u = ks1(upper);
        end
        if theta1>theta2
            bdu = [-pi,kesai_u];
        else
            bdu = [kesai_u, pi];
        end
        if lower<=l 
            kesai_l = kesai_s;
        else
            kesai_l = ks1(lower);  
        end
        if theta1<theta2
            bdl = [-pi,kesai_l];
        else
            bdl = [kesai_l, pi];
        end
        if (kesai_u - kesai_l)*(theta1-theta2)>0
            bd = bd_intersection(bdu, bdl);
        else
            bd = cat(2, bdl, bdu);
        end
         bounds = bd_intersection(bounds, bd);
    else
         if upper<=u && lower>=l
            bounds = [];
            return;
         end
        if upper<=u && upper>=l
            kesai_u = kesai_s;
        else
            kesai_u = ks1(upper);
        end
        if theta1>theta2
            bdu=[kesai_u, pi];
        else
            bdu=[-pi, kesai_u];
        end
        if lower>=l && lower<=u
            kesai_l = kesai_s;  
        else
            kesai_l = ks1(lower);
        end
        if theta1>theta2
            bdl=[-pi, kesai_l];
        else
            bdl=[kesai_l, pi];
        end
        if (kesai_u - kesai_l)*(theta1-theta2)<0
            bd = bd_intersection(bdu, bdl);
        else
            bd = cat(2, bdl, bdu);
        end
         bounds = bd_intersection(bounds, bd);   
    end
elseif delta > 0
    if bt == ct
        kesai1 = 0;
        kesai2 = pi;
    else
        kesai1 = 2*atan((at + sqrt(delta)) / (bt - ct));
        kesai2 = 2*atan((at - sqrt(delta)) / (bt - ct));
    end
    %       kesai1 = asin(-ct/(sqrt(at^2+bt^2))) - atan2(bt, at);
    %       kesai2 = pi - asin(-ct/(sqrt(at^2+bt^2))) - atan2(bt, at);
    if at*cos(kesai1)-bt*sin(kesai1) < 0 % local maxmia
        tem = kesai2;
        kesai2 = kesai1;
        kesai1 = tem;
    end
    % now kesai1: local minima, kesai2: local maxmia
    y1 = y(kesai1);
    y2 = y(kesai2);
    if abs(y1-y2) < 1e-6 % constant
        if upper >= y1 && lower <= y1
            bounds = [-pi, pi];
            return;
        else
            bounds=[];
            return;
        end
    end
    if y1 < y2 % no jump
        bounds = [-pi, pi];
        if lower >= y2
            bounds = [];
            return;
        elseif lower < y2 && lower > y1
            kesai = k2(lower);
            if d(kesai(1)) < 0
                bounds = [-pi, kesai(1), kesai(2), pi];
            else
                bounds = [kesai(1), kesai(2)];
            end
        end
        if upper <= y1
            bounds = [];
            return;
        elseif upper < y2 && upper > y1
            kesai = k2(upper);
            if d(kesai(1)) < 0
                bounds = bd_intersection(bounds, [kesai(1), kesai(2)]);
            else
                bounds = bd_intersection(bounds,[-pi, kesai(1), kesai(2), pi]);
            end
        end
    else % with jump
        bounds =[];
        if y1 < upper
            kesai = k2(upper);
            if (kesai(1)-kesai1)*(kesai(2)-kesai1) > 0
                bounds = [-pi, kesai(1), kesai(2), pi];
            else
                bounds = [kesai(1), kesai(2)];
            end
        end
        if y2 > lower
            kesai = k2(lower);
            if (kesai(1)-kesai2)*(kesai(2)-kesai2) > 0
                bounds = cat(2, bounds,[-pi, kesai(1), kesai(2), pi]);
            else
                bounds = cat(2, bounds,[kesai(1), kesai(2)]);
            end
        end
    end
else % delta < 0
    [kesai1, kesai2] = reorder(k1(lower), k1(upper));
    theta = y((kesai1 + kesai2) / 2);
    if theta > lower && theta < upper % no jump
        bounds = [kesai1, kesai2];
    else % with jump
        bounds = [-pi, kesai1, kesai2, pi];
    end
end
end

function kesais = anlysis_kesai(As, Bs, Cs, Aw, Bw, Cw, cfg, lowers, uppers, eps0)
[kesais2, kesai_s1] = kesai_range_hj('joint 2',  As(3,3), Bs(3,3), Cs(3,3), cfg(1),lowers(2), uppers(2),eps0);
if isempty(kesais2)
    kesais=[];
    return;
end
kesais1 = kesai_range_pj('joint 1', As(2,3), Bs(2,3), Cs(2,3), As(1,3), Bs(1,3), Cs(1,3), cfg(1),lowers(1), uppers(1),kesai_s1);
kesais3 = kesai_range_pj('joint 3', As(3,2), Bs(3,2), Cs(3,2), -As(3,1), -Bs(3,1), -Cs(3,1), cfg(1),lowers(3), uppers(3),kesai_s1);
kesais = bd_intersection(kesais1,kesais2);
kesais = bd_intersection(kesais,kesais3);
if isempty(kesais)
    if ~isempty(kesai_s1)
        kesais=[kesai_s1,kesai_s1];
    end
    return;
end
[kesais6, kesai_s2] = kesai_range_hj('joint 6',  Aw(3,3), Bw(3,3), Cw(3,3), cfg(3),lowers(6), uppers(6),eps0);
if isempty(kesais6)
    kesais=[];
    return;
end
kesais5 = kesai_range_pj('joint 5',  Aw(2,3), Bw(2,3), Cw(2,3), Aw(1,3), Bw(1,3), Cw(1,3), cfg(3),lowers(5), uppers(5),kesai_s2);
kesais7 = kesai_range_pj('joint 7',  Aw(3,2), Bw(3,2), Cw(3,2), -Aw(3,1), -Bw(3,1), -Cw(3,1), cfg(3),lowers(7), uppers(7),kesai_s2);
kesais567 = bd_intersection(kesais5,kesais6);
kesais567 = bd_intersection(kesais567,kesais7);
kesais = bd_intersection(kesais, kesais567);
if isempty(kesais) && ~isempty(kesai_s2)
    kesais=[kesai_s2,kesai_s2];
elseif isempty(kesais) && ~isempty(kesai_s1)
    kesais=[kesai_s1,kesai_s1];
end
end
