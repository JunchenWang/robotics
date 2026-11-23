function save_stereo_opencv_yaml(stereoParams, fileName)

    [K1, D1, K2, D2, R, T] = stereoParametersToOpenCV(stereoParams);

    K1 = double(K1);
    K2 = double(K2);

    D1 = reshape(double(D1), 1, []);
    D2 = reshape(double(D2), 1, []);

    R = double(R);
    T = reshape(double(T), 3, 1);

    fid = fopen(fileName, 'w');

    if fid < 0
        error("Cannot open file: %s", fileName);
    end

    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '%%YAML:1.0\n');
    fprintf(fid, '---\n');

    % fprintf(fid, 'image_width: %d\n', imageSize(2));
    % fprintf(fid, 'image_height: %d\n', imageSize(1));

    write_matrix(fid, 'K1', K1);
    write_matrix(fid, 'D1', D1);

    write_matrix(fid, 'K2', K2);
    write_matrix(fid, 'D2', D2);

    write_matrix(fid, 'R', R);
    write_matrix(fid, 'T', T);
end


function write_matrix(fid, name, matrix)

    [rows, cols] = size(matrix);

    fprintf(fid, '%s: !!opencv-matrix\n', name);
    fprintf(fid, '   rows: %d\n', rows);
    fprintf(fid, '   cols: %d\n', cols);
    fprintf(fid, '   dt: d\n');
    fprintf(fid, '   data: [ ');

    % OpenCV YAML 按行排列数据
    values = reshape(matrix.', 1, []);

    for i = 1:numel(values)
        if i ~= numel(values)
            fprintf(fid, '%.17g, ', values(i));
        else
            fprintf(fid, '%.17g', values(i));
        end
    end

    fprintf(fid, ' ]\n');
end