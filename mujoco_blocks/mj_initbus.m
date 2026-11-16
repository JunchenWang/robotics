function [a,b,c,d,e] = mj_initbus(xmlPath, cameraOutputEnabled)
% Copyright 2022-2023 The MathWorks, Inc.
% Lets just be on safe side and run mj_initbus_mex in a separate
% process. It calls glfw functions which work best in a main thread of
% a separate process.
persistent mh
if nargin < 2
    cameraOutputEnabled = true;
end
if ~(isa(mh,'matlab.mex.MexHost') && isvalid(mh))
    mh = mexhost;
end
[a,b,c,d,e] = feval(mh, 'mj_initbus_mex', xmlPath, logical(cameraOutputEnabled));
end
