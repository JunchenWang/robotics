function T = make_tform_qt(qt)
if isrow(qt)
    qt = qt';
end
T = [RMByQuaternion(qt(1:4)), qt(5:7); 0 0 0 1];