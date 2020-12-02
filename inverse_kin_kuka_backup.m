function [angles, bounds] = inverse_kin_kuka(R, t, cfg)
% inverse_kin_kuka kuka med的运动学逆解,自动选择kesai
% cfg signs of axis 2,4,6
lower = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
upper = -lower;
eps1 = 1e-8;
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

bounds = anlysis_kesai(As, Bs, Cs, Aw, Bw, Cw, cfg, lower, upper);
if isempty(bounds)
    return;
end
if length(bounds) == 2 || length(bounds) == 4 && bounds(2) - bounds(1) > bounds(4) - bounds(3)
    kesai = (bounds(1) + bounds(2)) / 2;
else
    kesai = (bounds(3) + bounds(4)) / 2;
end
if imag(kesai) ~= 0
    error('???');
end
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
angles = [theta1, theta2, theta3, theta4, theta5, theta6, theta7];
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
    x = linspace(-pi-0.01, pi+0.01, 1000);
%     x = x(2:end);
    figure;
    plot(x, y(x));
    title(str);
    grid on;
    delta = a^2+b^2;
    if abs(delta) < 1e-8 % a=0,b=0 theta is appx. constant wrt. kesasi
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
                if b*sin(ks(1)) - a*cos(ks(1)) > 0
                    bounds = bd_intersection(bounds, [-pi, ks(1), ks(2), pi]);
                else
                     bounds = bd_intersection(bounds, [ks(1), ks(2)]);
                end
%                 bounds(1) = ks(1);
%                 bounds(4) = ks(2);
            elseif cfg < 0 && theta2 < lower
                ks = k(lower);
                if b*sin(ks(1)) - a*cos(ks(1)) > 0
                    bounds = bd_intersection(bounds, [-pi, ks(1), ks(2), pi]);
                else
                     bounds = bd_intersection(bounds, [ks(1), ks(2)]);
                end
%                 bounds(1) = ks(1);
%                 bounds(4) = ks(2);
            end
        else
            d2 = (b*cos(kesai)+a*sin(kesai))/sin(theta);
            if d2 < 0 %maxima
               if cfg > 0 && theta2 >= upper || cfg < 0 && theta <= lower
                   bounds = []; 
               elseif cfg > 0 && theta > upper
                   ks = k(upper);
                   if b*sin(ks(1)) - a*cos(ks(1)) > 0
                       bounds= [-pi, ks(1), ks(2), pi];
                   else
                       bounds= [ks(1), ks(2)];
                   end
               elseif cfg < 0 && theta2 < lower
                   ks = k(lower);
                   if b*sin(ks(1)) - a*cos(ks(1)) > 0
                       bounds= [-pi, ks(1), ks(2), pi];
                   else
                       bounds= [ks(1), ks(2)];
                   end
               else
                   bounds = [-pi, pi];
               end 
            else % minima
               if cfg > 0 && theta >= upper || cfg < 0 && theta2 <= lower
                   bounds = [];  
               elseif cfg > 0 && theta2 > upper
                   ks = k(upper);
                   if b*sin(ks(1)) - a*cos(ks(1)) > 0
                       bounds= [-pi, ks(1), ks(2), pi];
                   else
                       bounds= [ks(1), ks(2)];
                   end
               elseif cfg < 0 && theta < lower
                   ks = k(lower);
                    if b*sin(ks(1)) - a*cos(ks(1)) > 0
                       bounds= [-pi, ks(1), ks(2), pi];
                   else
                       bounds= [ks(1), ks(2)];
                   end
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
    delta2 = bp^2 - 4*ap*cp;
    for kesai = [2*atan((-bp-sqrt(delta2)) / (2*ap)), 2*atan((-bp+sqrt(delta2)) / (2*ap))]
        if abs(kesai2theta_pj(kesai, an, bn, cn, ad, bd, cd, cfg) - theta) < 1e-2
            return;
        end
    end
    error('no solution');
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
    y = @(kesai) kesai2theta_pj(kesai, an, bn, cn, ad, bd, cd, cfg);
    k2 = @(theta) theta2kesai_pj_2(theta, an, bn, cn, ad, bd, cd);
    k1 = @(theta) theta2kesai_pj_1(theta, an, bn, cn, ad, bd, cd, cfg);
    x = linspace(-pi-0.01, pi+0.01, 1000);
%     x = x(2:end);
    figure;
    plot(x, y(x));
    title(str);
    grid on;
    % theta = atan2(u, v)
    at = cn*bd-bn*cd; %原论文有错，不带cfg
    bt = an*cd-cn*ad;
    ct = an*bd-bn*ad;
    d = @(kesai) at*sin(kesai)+bt*cos(kesai)+ct;
    if at==0 && bt==0 && ct==0 % theta is constant wrt. kesai
        theta = y(0);
        if upper >= theta && lower <= theta
            bounds = [-pi, pi];
            return;
        else
            bounds=[];
            return;
        end
    end
    delta = (at^2 + bt^2 - ct^2);
    delta_n = delta /(at^2 + bt^2 + ct^2);
%     disp(delta/(at^2 + bt^2 + ct^2));
    if abs(delta_n) <= 1e-3
        tol = 0.1; % 5.7 degree
        if abs(bt-ct)/abs(bt) > 1e-2
            kesai_s = 2*atan(at/(bt-ct)); % singular arm angle
        else
%             error('check it 1!');
              kesai_s = 2*atan(-(bt+ct)/(2*at));
              if kesai_s > 0
                  kesai2 = kesai_s - tol;
                  kesai1 = kesai_s - 2*pi + tol;
              else
                  kesai1 = kesai_s + tol;
                  kesai2 = kesai_s + 2*pi - tol;
              end
              theta1 = y(kesai1);
              theta2 = y(kesai2);
              bounds = [kesai1, kesai2];
              if (theta2-theta1)*d(kesai1) > 0 % no jump
                  u = max(theta1, theta2);
                  l = min(theta1, theta2);
                  if upper <= l || lower >= u
                      bounds = [];
                      return;
                  end
                  if upper < u
                      for kesai = k(upper)
                          if abs(y(kesai) - upper) < 1e-3
                              break;
                          end
                      end
                      if d(kesai) > 0
                          bounds = bd_intersection(bounds, [kesai1, kesai]);
                      else
                          bounds = bd_intersection(bounds, [kesai, kesai2]);
                      end
                  end
                  if lower > l
                      for kesai = k(lower)
                          if abs(y(kesai) - lower) < 1e-3
                              break;
                          end
                      end
                      if d(kesai) > 0
                          bounds = bd_intersection(bounds, [kesai, kesai2]);
                      else
                          bounds = bd_intersection(bounds, [kesai1, kesai]);
                      end
                  end
              else % with jump
                   for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-2
                            break;
                        end
                   end
                   kesai1 = kesai;
                   for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-2
                            break;
                        end
                   end
                   kesai2 = kesai;
                   if kesai1 > kesai2
                       tem = kesai1;
                       kesai1 = kesai2;
                       kesai2 = tem;
                   end
                   bounds = bd_intersection(bounds, [-pi,kesai1, kesai2, pi]);
              end
              
              return;
        end
        kesai1 = kesai_s - tol;
        kesai2 = kesai_s + tol;
        bounds(1) = -pi;
        bounds(2) = kesai1;
        bounds(3) = kesai2;
        bounds(4) = pi;
        theta1 = y(kesai1);
        d1 = d(kesai1);
        theta2 = y(kesai2);
%         d2 = at*sin(kesai2)+bt*cos(kesai2)+ct;  
        if theta1 > theta2 % jump -pi
            if d1 > 0 % 1
                if upper <= theta2 || lower >= theta1
                    bounds = [];
                    return;
                end
                if theta1 > upper
                    for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-3
                            break;
                        end
                    end
                    if kesai < kesai1
                        bounds(2) = kesai;
                    else
                        bounds = [kesai2, kesai];
                    end
                end
                if theta2 < lower
                    for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-3
                            break;
                        end
                    end
                    if kesai > kesai2
                        bounds = bd_intersection(bounds, [-pi, kesai1, kesai, pi]);
                    else
                        bounds = bd_intersection(bounds, [-pi, kesai, kesai2, pi]);
                    end
                end
            else % 2
                flag1 = 0; flag2 = 0;
                if theta1 < upper
                    for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-3
                            break;
                        end
                    end
                    if kesai < kesai1
                        bounds(1) = kesai;
                    end
                else
                    flag1 = 1;
                end
                if theta2 > lower
                    for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-3
                            break;
                        end
                    end
                    if kesai > kesai2
                        bounds(4) = kesai;
                    end
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
           if d1 < 0 % 3
                if upper <= theta1 || lower >= theta2
                    bounds = [];
                    return;
                end
                if theta1 < lower
                    for kesai = k(lower)
                        if abs(y(kesai) - lower) < 1e-3
                            break;
                        end
                    end
                    if kesai < kesai1
                        bounds(2) = kesai;
                    else
                        bounds = [kesai2, kesai];
                    end
                end
                if theta2 > upper
                    for kesai = k(upper)
                        if abs(y(kesai) - upper) < 1e-3
                            break;
                        end
                    end
                    if kesai > kesai2
                        bd_intersection(bounds, [-pi, kesai1, kesai, pi]);
                    else
                        bd_intersection(bounds, [kesai, kesai1]);
                    end
                end
            else % 4
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
    elseif delta_n > 0
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
            if lower >= y2
                bounds = [];
                return;
            elseif lower < y2 && lower > y1
                kesai = k(lower);
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
                kesai = k(upper);
                if d(kesai(1)) < 0
                    bounds = bd_intersection(bounds, [kesai(1), kesai(2)]);
                else
                    bounds = bd_intersection(bounds,[-pi, kesai(1), kesai(2), pi]);
                end
            end
        else % with jump
            bounds =[];
            if kesai1 < kesai2
                if y1 < upper
                    kesai = k(upper);
                    if kesai(1) > kesai1
                        bounds = [-pi, kesai(1)];
                    else
                        bounds = [kesai(1), kesai(2)];
                    end
                end
                if y2 > lower
                    kesai = k(lower);
                    if kesai(2) < kesai2
                        bounds(3:4) = [kesai(2), pi];
                    else
                        bounds(3:4) = [kesai(1), kesai(2)];
                    end
                end
            else
                if y1 < upper
                    kesai = k(upper);
                    if kesai(2) < kesai1
                        bounds = [kesai(2), pi];
                    else
                        bounds = [kesai(1), kesai(2)];
                    end
                end
                if y2 > lower
                    kesai = k(lower);
                    if kesai(1) > kesai2
                        bounds(3:4) = [-pi, kesai(1)];
                    else
                        bounds(3:4) = [kesai(1), kesai(2)];
                    end
                end
            end
        end
    else % delta_n < 0
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
kesais3 = kesai_range_pj('joint 3', As(3,2), Bs(3,2), Cs(3,2), -As(3,1), -Bs(3,1), -Cs(3,1), cfg(1),lower(3), upper(3));% A3
kesais5 = kesai_range_pj('joint 5',  Aw(2,3), Bw(2,3), Cw(2,3), Aw(1,3), Bw(1,3), Cw(1,3), cfg(3),lower(5), upper(5));% A1
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
