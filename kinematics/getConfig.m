function cfg = getConfig(angles)
cfg = ones(1,3);
if angles(2) < 0
    cfg(1) = -1;
end
if angles(4) < 0
    cfg(2) = -1;
end
if angles(6) < 0
    cfg(3) = -1;
end