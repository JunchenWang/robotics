function [mass, r, offset] = force_calibration(forceFile, poseFile)
fid = fopen(forceFile,'r'); % fx fy fz mx my mz
force = fscanf(fid, '%f', [6, inf]);
fclose(fid);
fid = fopen(poseFile,'r'); % x y z rx ry rz
pose = fscanf(fid, '%f', [6, inf]);
fclose(fid);
data = [pose', force'];
[mass, r, offset] = force_sensor_calibration(data);