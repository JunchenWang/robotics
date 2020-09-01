function qd = ConcatenateDualQuaternions(qd1, qd2, varargin)
if nargin == 2
qd = [ConcatenateQuaternions(qd1(1 : 4), qd2(1 : 4));
      ConcatenateQuaternions(qd1(1 : 4), qd2(5 : 8)) + ConcatenateQuaternions(qd1(5 : 8), qd2(1 : 4))];
else
    qd = ConcatenateDualQuaternions(qd1, qd2);
    n = length(varargin);
    for i = 1 : n
        qd = ConcatenateDualQuaternions(qd, varargin{i});
    end
end
