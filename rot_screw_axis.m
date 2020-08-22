function sa = rot_screw_axis(p, dir)
% v is nx6 each row represents a twist
sa = screw_axis(p, dir, 0);