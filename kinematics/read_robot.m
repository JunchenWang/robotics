function [robot, joint_names] = read_robot(source, varargin)
%READ_ROBOT Read a robot chain from a URDF or MuJoCo MJCF file.
%   ROBOT = READ_ROBOT(SOURCE, FLANGE_LINK, BASE_LINK) selects the reader
%   from the XML root. The optional endpoints use the selected reader's
%   defaults when omitted. SOURCE may also be URDF XML text.
%   [ROBOT, JOINT_NAMES] = READ_ROBOT(...) also returns the movable joints
%   in base-to-flange order.

narginchk(1, 3);
validateattributes(source, {'char', 'string'}, {'scalartext', 'nonempty'}, 'read_robot', 'source');
source = char(source);
isXmlText = startsWith(strtrim(source), '<');
parser = matlab.io.xml.dom.Parser;
% Match read_urdf's handling of extension tags with undeclared prefixes.
parser.Configuration.Namespaces = false;
if isXmlText
    document = parseString(parser, source);
else
    filename = source;
    if ~isfile(filename)
        filename = which(filename);
    end
    if isempty(filename)
        error('read_robot:FileNotFound', 'Robot file not found: %s', source);
    end
    source = filename;
    document = parseFile(parser, source);
end

root = document.getDocumentElement();
tag = root.getTagName();
switch tag
    case 'robot'
        reader = @read_urdf;
    case 'mujoco'
        if isXmlText
            error('read_robot:MJCFFileRequired', 'Specify a file path for an MJCF model.');
        end
        reader = @read_mujoco;
    otherwise
        error('read_robot:UnsupportedRoot', 'Unsupported XML root <%s>. Expected <robot> or <mujoco>.', tag);
end
[robot, joint_names] = reader(source, varargin{:});
end
