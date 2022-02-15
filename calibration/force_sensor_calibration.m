function [mass, r, offset] = force_sensor_calibration(data)
n = size(data,1);
m = n*(n-1) / 2;
M = zeros(1, m);
offsets = zeros(6,m);
A = zeros(3 * m, 3);
b = zeros(3 * m, 1);
cnt = 1;
for i = 1 : n - 1
    for j = i + 1 : n
        Ri = RotationByAxisAngleRep(data(i,4:6));
        Rj = RotationByAxisAngleRep(data(j,4:6));
        Fi = data(i,7:9);
        Fj = data(j,7:9);
        M(cnt) = norm(Fi - Fj) / norm((Ri - Rj)'*[0,0,-1]');
        cnt = cnt + 1;
    end
end
plot(M);
mass = mean(M);
cnt = 1;
for i = 1 : n - 1
    for j = i + 1 : n
        Ri = RotationByAxisAngleRep(data(i,4:6));
        Rj = RotationByAxisAngleRep(data(j,4:6));
        Mi = data(i,10:12);
        Mj = data(j,10:12);
        b(3 * cnt - 2 : 3 * cnt) = -(Mi - Mj)';
        A(3 * cnt - 2 : 3 * cnt, :) = so_w((Ri - Rj)' * [0,0,-1]') * mass;
        cnt = cnt + 1;
    end
end
r = A \ b;
if (norm(b - A * r) / norm(b) > 0.02)
    warning('data may be not good');
end
cnt = 1;
for i = 1 : n - 1
    for j = i + 1 : n
        Ri = RotationByAxisAngleRep(data(i,4:6));
        Rj = RotationByAxisAngleRep(data(j,4:6));
        Fi = data(i,7:9);
        Fj = data(j,7:9);
        Mi = data(i,10:12);
        Mj = data(j,10:12);
        offsets(1:3, cnt) = (Fi' + Fj' - mass * (Ri + Rj)' * [0, 0, -1]') / 2;
        offsets(4:6, cnt) = (Mi' + Mj' - mass * cross(r, (Ri + Rj)' * [0, 0, -1]')) / 2;
        cnt = cnt + 1;
    end
end
hold on;
plot(offsets');
offset = mean(offsets,2);
fid = fopen('sensor_calib.txt', 'w');
fprintf(fid, '%f ', [mass, r', offset']);
fclose(fid);