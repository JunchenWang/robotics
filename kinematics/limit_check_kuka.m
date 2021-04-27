function flag = limit_check_kuka(angles)
lowers = [-170, -120, -170, -120 ,-170, -120, -175] / 180 * pi;
uppers = -lowers;
flag = any(~(angles >= lowers & angles <= uppers));