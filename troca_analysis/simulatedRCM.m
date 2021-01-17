function data = simulatedRCM
rcm = [20, 20, 20];
N = 10;
data = zeros(N, 6);
for i = 1 : N
    dir = rand(1,3);
    dir = dir / norm(dir);
    ang = 0.2 + rand();
    pitch = (100 * rand() - 50) / ang;
    sm = [rcm, dir, pitch, ang];
    twist = sm2twist(sm);
    plot_twist(twist);
    T = exp_twist(twist);
    data(i,:) = RVFrame2EulZYXFrame(Frame_T_Converter(T));
end