function hand_eye_calib
fidrobot = fopen('D:\MICR\MICSys\temp\hand-eye\robot_pos.txt', 'r');
robot = fscanf(fidrobot, '%f', [6, inf]);
fidmarker = fopen('D:\MICR\MICSys\temp\hand-eye\marker_pos.txt', 'r');
marker = fscanf(fidmarker, '%f', [6, inf]);

[fX, frame_err] = HandEyeCalibration(robot, marker, 1);
disp(fX);
disp(frame_err);