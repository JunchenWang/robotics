%Save stereo calibration results
function saveStereoCalibResults(stereoParams)
xdoc=com.mathworks.xml.XMLUtils.createDocument('opencv_storage');
xroot=xdoc.getDocumentElement;

%========定义矩阵，写入xml
%%左相机
%内参
type=xdoc.createElement('cameraMatrix_left');
xroot.appendChild(type);
type.setAttribute('type_id','opencv-matrix')

cameraMatrix_left=stereoParams.CameraParameters1.IntrinsicMatrix;
cameraMatrix_left(3,1)=cameraMatrix_left(3,1)-1;
cameraMatrix_left(3,2)=cameraMatrix_left(3,2)-1;
[m1,n1]=size(cameraMatrix_left);

rows=xdoc.createElement('rows');
rows.appendChild(xdoc.createTextNode(sprintf('%d',m1)));
type.appendChild(rows);

cols=xdoc.createElement('cols');
cols.appendChild(xdoc.createTextNode(sprintf('%d',n1)));
type.appendChild(cols);

dt=xdoc.createElement('dt');
dt.appendChild(xdoc.createTextNode(sprintf('%s','d')));
type.appendChild(dt);

data=xdoc.createElement('data');
data.appendChild(xdoc.createTextNode(sprintf('%f ',cameraMatrix_left)));
type.appendChild(data);
%畸变
type=xdoc.createElement('distCoeffs_left');
xroot.appendChild(type);
type.setAttribute('type_id','opencv-matrix')

distCoeffs_left=[stereoParams.CameraParameters1.RadialDistortion(1),stereoParams.CameraParameters1.RadialDistortion(2),stereoParams.CameraParameters1.TangentialDistortion(1),stereoParams.CameraParameters1.TangentialDistortion(2)];

[m1,n1]=size(distCoeffs_left);

rows=xdoc.createElement('rows');
rows.appendChild(xdoc.createTextNode(sprintf('%d',m1)));
type.appendChild(rows);

cols=xdoc.createElement('cols');
cols.appendChild(xdoc.createTextNode(sprintf('%d',n1)));
type.appendChild(cols);

dt=xdoc.createElement('dt');
dt.appendChild(xdoc.createTextNode(sprintf('%s','d')));
type.appendChild(dt);

data=xdoc.createElement('data');
data.appendChild(xdoc.createTextNode(sprintf('%f ',distCoeffs_left)));
type.appendChild(data);
%%右相机
%内参
type=xdoc.createElement('cameraMatrix_right');
xroot.appendChild(type);
type.setAttribute('type_id','opencv-matrix')

cameraMatrix_right=stereoParams.CameraParameters2.IntrinsicMatrix;
cameraMatrix_right(3,1)=cameraMatrix_right(3,1)-1;
cameraMatrix_right(3,2)=cameraMatrix_right(3,2)-1;
[m1,n1]=size(cameraMatrix_right);

rows=xdoc.createElement('rows');
rows.appendChild(xdoc.createTextNode(sprintf('%d',m1)));
type.appendChild(rows);

cols=xdoc.createElement('cols');
cols.appendChild(xdoc.createTextNode(sprintf('%d',n1)));
type.appendChild(cols);

dt=xdoc.createElement('dt');
dt.appendChild(xdoc.createTextNode(sprintf('%s','d')));
type.appendChild(dt);

data=xdoc.createElement('data');
data.appendChild(xdoc.createTextNode(sprintf('%f ',cameraMatrix_right)));
type.appendChild(data);
%畸变
type=xdoc.createElement('distCoeffs_right');
xroot.appendChild(type);
type.setAttribute('type_id','opencv-matrix')

distCoeffs_right=[stereoParams.CameraParameters2.RadialDistortion(1),stereoParams.CameraParameters2.RadialDistortion(2),stereoParams.CameraParameters2.TangentialDistortion(1),stereoParams.CameraParameters2.TangentialDistortion(2)];

[m1,n1]=size(distCoeffs_right);

rows=xdoc.createElement('rows');
rows.appendChild(xdoc.createTextNode(sprintf('%d',m1)));
type.appendChild(rows);

cols=xdoc.createElement('cols');
cols.appendChild(xdoc.createTextNode(sprintf('%d',n1)));
type.appendChild(cols);

dt=xdoc.createElement('dt');
dt.appendChild(xdoc.createTextNode(sprintf('%s','d')));
type.appendChild(dt);

data=xdoc.createElement('data');
data.appendChild(xdoc.createTextNode(sprintf('%f ',distCoeffs_right)));
type.appendChild(data);
%%旋转矩阵，Matlab和opencv的存储方式不同，
%但是由于二者左乘右乘的方式不同，因此直接保存就可以
type=xdoc.createElement('R');
xroot.appendChild(type);
type.setAttribute('type_id','opencv-matrix')

R=stereoParams.RotationOfCamera2;

[m1,n1]=size(R);

rows=xdoc.createElement('rows');
rows.appendChild(xdoc.createTextNode(sprintf('%d',m1)));
type.appendChild(rows);

cols=xdoc.createElement('cols');
cols.appendChild(xdoc.createTextNode(sprintf('%d',n1)));
type.appendChild(cols);

dt=xdoc.createElement('dt');
dt.appendChild(xdoc.createTextNode(sprintf('%s','d')));
type.appendChild(dt);

data=xdoc.createElement('data');
data.appendChild(xdoc.createTextNode(sprintf('%f ',R)));
type.appendChild(data);

%%平移矩阵
type=xdoc.createElement('T');
xroot.appendChild(type);
type.setAttribute('type_id','opencv-matrix')

T=stereoParams.TranslationOfCamera2;

[m1,n1]=size(T);

rows=xdoc.createElement('rows');
rows.appendChild(xdoc.createTextNode(sprintf('%d',n1)));
type.appendChild(rows);

cols=xdoc.createElement('cols');
cols.appendChild(xdoc.createTextNode(sprintf('%d',m1)));
type.appendChild(cols);

dt=xdoc.createElement('dt');
dt.appendChild(xdoc.createTextNode(sprintf('%s','d')));
type.appendChild(dt);

data=xdoc.createElement('data');
data.appendChild(xdoc.createTextNode(sprintf('%f ',T)));
type.appendChild(data);

str=strcat('stereo_calibration9_7','.xml');
xmlwrite(str,xdoc);
end