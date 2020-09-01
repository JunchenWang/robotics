function tf = Frame_T_Converter(ft)
if isvector(ft)
    tf = eye(4);
    tf(1 : 3, 1 : 3) = RotationByAxisAngleRep(ft(1 : 3));
    tf(1 : 3, 4) = ft(4 : 6); 
else
    tf = zeros(6, 1);
    tf(1 : 3) = AngleAxisFromRotation(ft(1 : 3, 1 : 3));
    tf(4 : 6) = ft(1 : 3, 4);
end