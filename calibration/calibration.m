% function [mass, r, offset] = calibration(forceFile, poseFile)
fid = fopen('force.txt','r');
force = fscanf(fid, '%f', [6, inf]);
fclose(fid);
fid = fopen('pose.txt','r');
pose = fscanf(fid, '%f', [6, inf]);
fclose(fid);
data = [pose', force'];
[mass, r, offset] = force_sensor_calibration(data);