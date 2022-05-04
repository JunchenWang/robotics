
function impedence_control_simulation

robot = read_dynamics_file('F:\robotics\urdf\iiwa7\dynamics.txt');
% robot = UR5e;
u = udpport("byte");
ptp([-40, 70, 0, -80, 0, -60, 0]/180*pi);
% ptp([0, -60, 80, -100, -90, 0]/180*pi);

lineTo2(robot, [-500,0,0]); % axang2rotm([0,1,0,pi/2]));

    function F = Wrench(t)
        if t < 3 && t > 1
            F = [12, 0, 0, 0, 0, 10]';
        else
            F = zeros(6,1);
        end
    end

    function lineTo2(robot, t, R)
        if nargin < 3
            R = eye(3);
        end
        if isrow(t)
            t = t';
        end
        Tcp = eye(4);
        Tsensor = eye(4);
        start = queryJoints;
        pre_q = start;
        q = pre_q;
        kesai = cal_kuka_kesai(start);
        Ts = forward_kin_kuka(start);
        Ts(1:3,4) = Ts(1:3,4) / 1000;
        Te = Ts*[R,t / 1000;0 0 0 1];
        T = 4;
        Freq = 500;
        dt = 1 / Freq;
        %     r = rateControl(Freq);
        numSamples = round(T * Freq) + 1;
        data = zeros(robot.dof, numSamples);
        [s,sd,sdd] = trapveltraj([0, 1],numSamples, 'EndTime', T);
        tSamples = linspace(0,T,numSamples);
        [tforms,vel,~] = transformtraj(Ts,Te,[0 T],tSamples, 'TimeScaling', [s;sd;sdd]);
        m = 1;
        k1 = 400;
        k2 = 40;
        b1 = 2*sqrt(k1*m);
        b2 = 2*sqrt(k2*m);
        M = cat(3, m*eye(3), m*eye(3));
        B = cat(3, b1 * eye(3), b2 * eye(3));
        K = cat(3, k1 * eye(3), k2 * eye(3));
        re2 = zeros(3,1);
        pe2 = zeros(3,1);
        red2 = zeros(3,1);
        ped2 = zeros(3,1);
                re = zeros(3,1);
        pe = zeros(3,1);
        red = zeros(3,1);
        ped = zeros(3,1);
        data_re = zeros(6, numSamples);
        data_pe = zeros(6, numSamples);
        for i = 1 : numSamples
            Td = tforms(:,:,i);
%             Td = Ts;
%             Vd = zeros(6,1);
            Vd = vel(:, i);
            qd = (q - pre_q) / dt;
            data(:,i) = qd';
            y = [q,qd]';
            f = Wrench(tSamples(i));
            pedd2 = M(:,:,1) \ (f(4:6) - B(:,:,1) * ped2 - K(:,:,1) * pe2);
            redd2 = M(:,:,2) \ (f(1:3) - B(:,:,2) * red2 - K(:,:,2) * re2);
            ped2 = ped2 + dt * pedd2;
            pe2 = pe2 + dt * ped2;
            red2 = red2 + dt * redd2;
            re2 = re2 + dt * red2;
% 
%               [re, pe] = impedence_control(robot, Td, Vd, M(:,:,1), B(:,:,1), K(:,:,1),...
%                                                  M(:,:,2), B(:,:,2), K(:,:,2), y, f, dt);
[T, re, pe, red, ped, redd, pedd] = impedance_control(robot, Tcp, Tsensor, Td, Vd, M(:,:,1), B(:,:,1), K(:,:,1),...
                                                   M(:,:,2), B(:,:,2), K(:,:,2), y, f, dt);
            data_pe(:,i) = [pe; 0;0;0];
            data_re(:,i) = [re; 0;0;0];
            Tbs = eye(4);
            Tbs(1:3,1:3) = Td(1:3,1:3)';
            Vd_b = adjoint_T(Tbs) * Vd * dt;
            Td = Td * exp_twist(Vd_b);
%             delta = Td - tforms(:,:,i+1);
%             disp(norm(delta));
            R = Td(1:3,1:3) * exp_w(-re);
            p = Td(1:3,4) - R * pe;
            X = [R, p; 0, 0, 0, 1];
            angles = inverse_kin_kuka_kesai_near(X, kesai, q);
%             angles = inverse_kin_general(robot, X, q, [1e-5,1e-5]);
            %         XX = forward_kin_general(robot ,angles);
            %         disp(norm(XX - X));
            if isempty(angles)
                plot(data');
                error('no solution');
            end
            setJoints(angles);
            pre_q = q;
            q = angles;
            %           [re, pe] = impedence_control(robot, Te, zeros(6,1), M, B, K, y, f, dt, re, pe);
            %         angles = inverse_kin_kuka_kesai_near(X, kesai, q);
            %         [angles, flag] = UR_inverse_kin_near(robot, X, pre_q);
            %         waitfor(r);
        end
        plot(data_pe');
        figure;
        plot(data_re');
    end

    function ptp(jts, vel)
        if nargin < 2
            vel = 1;
        end
        start = queryJoints;
        wayPoints = [start',jts'];
        Freq = 200;
        r = rateControl(Freq);
        T = max(abs(jts - start) / vel);
        numSamples = round(T * Freq) + 1;
        [q,~,~,~,~] = trapveltraj(wayPoints,numSamples);

        for i = 1 : numSamples
            setJoints(q(:,i));
            waitfor(r);
        end
    end

    function joints = queryJoints
        % ";" 表示查询关节角
        writeline(u,"robot;","127.0.0.1",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
    end

    function gohome
        ptp([0,0,0,0,0,0,0]);
    end

    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2), jt(3), jt(4), jt(5)...
            ,jt(6), jt(7));
        writeline(u,cmd,"127.0.0.1",7755);
    end
end

