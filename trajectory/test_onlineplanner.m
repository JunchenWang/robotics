function test_onlineplanner

fid = fopen('F:\MICR\launcher\data.txt');
data = fscanf(fid, '%f', [18, inf])';
t = 0.001 * (0 : size(data, 1) - 1);
plot(t, data, 'LineWidth',2);

