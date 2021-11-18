function [score, T] = imageMetricUnderTransform(original, distorted, thr)

ptsOriginal  = detectSURFFeatures(original, 'MetricThreshold', thr);
ptsDistorted = detectSURFFeatures(distorted, 'MetricThreshold', thr);

[featuresOriginal,  validPtsOriginal]  = extractFeatures(original,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(distorted, ptsDistorted);

indexPairs = matchFeatures(featuresOriginal, featuresDistorted);
matchedOriginal  = validPtsOriginal.Location(indexPairs(:,1),:);
matchedDistorted = validPtsDistorted.Location(indexPairs(:,2),:);

[tform, ~] = estimateGeometricTransform2D(...
    matchedDistorted, matchedOriginal, 'similarity');

T = tform.T;

outputView = imref2d(size(original));
recovered  = imwarp(distorted,tform,'OutputView',outputView);

score = sum((single(original) - single(recovered)).^2, 'all') / numel(original);
% imshowpair(original,recovered,'montage')


