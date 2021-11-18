original = im2gray(imread('2.bmp'));
distorted = im2gray(imread('1.bmp'));
imshowpair(original,distorted,'montage')


%% Step 3: Find Matching Features Between Images
% Detect features in both images.
ptsOriginal  = detectSURFFeatures(original, 'MetricThreshold', 100);
ptsDistorted = detectSURFFeatures(distorted, 'MetricThreshold', 100);

%%
% Extract feature descriptors.
[featuresOriginal,  validPtsOriginal]  = extractFeatures(original,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(distorted, ptsDistorted);

%%
% Match features by using their descriptors.
indexPairs = matchFeatures(featuresOriginal, featuresDistorted);

%%
% Retrieve locations of corresponding points for each image.
matchedOriginal  = validPtsOriginal(indexPairs(:,1));
matchedDistorted = validPtsDistorted(indexPairs(:,2));

%%
% Show putative point matches.
figure;
showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
title('Putatively matched points (including outliers)');

%% Step 4: Estimate Transformation
% Find a transformation corresponding to the matching point pairs using the
% statistically robust M-estimator SAmple Consensus (MSAC) algorithm, which
% is a variant of the RANSAC algorithm. It removes outliers while computing
% the transformation matrix. You may see varying results of the
% transformation computation because of the random sampling employed by the
% MSAC algorithm.
[tform, inlierIdx] = estimateGeometricTransform2D(...
    matchedDistorted, matchedOriginal, 'similarity');
inlierDistorted = matchedDistorted(inlierIdx, :);
inlierOriginal  = matchedOriginal(inlierIdx, :);
%%
% Display matching point pairs used in the computation of the
% transformation.
figure; 
showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted);
title('Matching points (inliers only)');
legend('ptsOriginal','ptsDistorted');

%% Step 5: Solve for Scale and Angle
% Use the geometric transform, tform, to recover the scale and angle.
% Since we computed the transformation from the distorted to the original
% image, we need to compute its inverse to recover the distortion.
%
%  Let sc = s*cos(theta)
%  Let ss = s*sin(theta)
%
%  Then, Tinv = [sc -ss  0;
%                ss  sc  0;
%                tx  ty  1]
%
%  where tx and ty are x and y translations, respectively.
%

%%
% Compute the inverse transformation matrix.
Tinv  = tform.invert.T;

ss = Tinv(2,1);
sc = Tinv(1,1);
scaleRecovered = sqrt(ss*ss + sc*sc)
thetaRecovered = atan2(ss,sc)*180/pi

%%
% The recovered values should match your scale and angle values selected in
% *Step 2: Resize and Rotate the Image*.

%% Step 6: Recover the Original Image
% Recover the original image by transforming the distorted image.
outputView = imref2d(size(original));
recovered  = imwarp(distorted,tform,'OutputView',outputView);

%%
% Compare |recovered| to |original| by looking at them side-by-side in a
% montage.
figure, imshowpair(original,recovered,'montage')


