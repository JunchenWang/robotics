function hand_eye_calib
fidrobot = fopen('robot_pos.txt', 'r');
robot = fscanf(fidrobot, '%f', [6, inf]);
fidmarker = fopen('marker_pos.txt', 'r');
marker = fscanf(fidmarker, '%f', [6, inf]);

[fX, frame_err] = HandEyeCalibration(robot, marker, 1);
disp(fX);
disp(frame_err);
fclose(fidrobot);
fclose(fidmarker);