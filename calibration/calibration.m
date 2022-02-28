function [mass, r, offset] = calibration(forceFile, poseFile)
fid = fopen(forceFile,'r');
force = fscanf(fid, '%f', [6, inf]);
fclose(fid);
fid = fopen(poseFile,'r');
pose = fscanf(fid, '%f', [6, inf]);
fclose(fid);
data = [pose', force'];
[mass, r, offset] = force_sensor_calibration(data);