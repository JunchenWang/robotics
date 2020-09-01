function q = ConcatenateQuaternions(q2, q1, varargin)
if nargin == 2
    q2 = q2(:);
    q1 = q1(:);
    a1 = q1(1);
    a2 = q2(1);
    v1 = q1(2 : 4);
    v2 = q2(2 : 4);
    q = [a1 * a2 - dot(v1, v2); cross(v2, v1) + a2 * v1 + a1 * v2];
else
    n = length(varargin);
    q = ConcatenateQuaternions(q2, q1);
    for i = 1 : n
        q = ConcatenateQuaternions(q, varargin{i});
    end 
end