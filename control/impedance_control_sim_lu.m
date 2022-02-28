function impedance_control_sim_lu

robot = read_dynamics_file('F:\robotics\urdf\iiwa7\dynamics.txt');
deg_f=7;

u = udpport("byte");
jtarget=[0, 70, 0, -80, 0, -60, 0]'/180*pi;
Td=forward_kin_kuka(jtarget);
Td(1:3,4)=Td(1:3,4)/1000;
jstart=jtarget;
ptp(jstart');
T=forward_kin_kuka(jstart);
T(1:3,4)=T(1:3,4)/1000;
   Rc=T(1:3,1:3);
tc=T(1:3,4);
Rd=Td(1:3,1:3);
td=Td(1:3,4);



q=jstart;
qd=zeros(deg_f,1);
vd=zeros(6,1);


Mr=100*eye(3,3);
Mt=100*eye(3,3);
Dr=10*eye(3,3);
Dt=10*eye(3,3);
Kr=300*eye(3,3);
Kt=300*eye(3,3);

J=jacobian_matrix(robot, q);
X_Re=logR(Rc'*Rd)';
X_Red=A(X_Re)\sym2vec(Rd'*Rc*(SkewMatrix(J(1:3,:)*qd)'*Rc'*Rd+Rc'*SkewMatrix(vd(1:3))*Rd));

X_te=Rc'*(td-tc);
X_ted=SkewMatrix(J(1:3,:)*qd)'*Rc'*(td-tc)+Rc'*vd(4:6)-J(4:6,:)*qd;
dt=0.002;
time =2;
n=2;
%trajectory
tInterval = [0 time];
tvec = 0:dt:time+dt;
% delta=[eye(3,3),[0;0;-0.2];0,0,0,1];
% Tfinal=T*delta;
Tfinal=T;
Tfinal(3,4)=Tfinal(3,4)+0.4;
[tfInterp, vd, ~] = transformtraj(T,Tfinal,tInterval,tvec);
% for i=1:size(vd,2)
%   [angles, ~] =UR_inverse_kin_near(robot, tfInterp(:,:,i), q);
%     setJoints(angles);
% end


for j=0:dt:time
  if j<5
    tao=[0;0;0];
    F=[0;0;0];
  else
      tao=[0;0;0];
     F=zeros(3,1);
  end
   

%    tspan = [0 dt];
%     yr0= [X_Re;X_Red];
%     [t,y] = ode45(@(t,y)impendance(t,y,Mr,Dr,Kr,tao), tspan, yr0);
%     X_Re=y(end,1:3)';
%     
%     yt0=[X_te;X_ted];
%     [t,y] = ode45(@(t,y)impendance(t,y,Mt,Dt,Kt,F), tspan, yt0);
%     X_te=y(end,1:3)';
  
  
  X_Redd=Mr\(tao-Kr*X_Re-Dr*X_Red);
  X_tedd=Mt\(F-Kt*X_te-Dt*X_ted);
    %%%更新X
  X_Re=X_Re+X_Red*dt+0.5*X_Redd*dt^2;
  X_te=X_te+X_ted*dt+0.5*X_tedd*dt^2;

  Td=tfInterp(:,:,n);
  Rd=Td(1:3,1:3);
  td=Td(1:3,4);
  Rc=Rd*exp_w(X_Re)';
  tc=td-Rc*X_te;
  T=[Rc,tc;0,0,0,1];
  kesai = cal_kuka_kesai(q);
  angles =inverse_kin_kuka_kesai_near(T, kesai, q);
  if size(angles,2)>1
      angles=angles';
  end
   if isempty(angles)
            error('no solution');
   end
   
    %angles=angles';

    if(norm(angles-q)>1)    
        disp("angles");
    end
    qd=(angles-q)/dt;
    q=angles;
    x(n-1,:)=q;
    setJoints(q');
    
    [J, ~] = jacobian_matrix(robot, q);
    X_Red=A(X_Re)\sym2vec(Rd'*Rc*(SkewMatrix(J(1:3,:)*qd)'*Rc'*Rd+Rc'*SkewMatrix(vd(1:3,n))*Rd));
    X_ted=SkewMatrix(J(1:3,:)*qd)'*Rc'*(td-tc)+Rc'*vd(4:6,n)-J(4:6,:)*qd;
    
    n=n+1;
end
plot(x)

    

    function setJoints(jt)
        cmd = sprintf('robot;%f;%f;%f;%f;%f;%f;%f;', jt(1), jt(2),jt(3),jt(4), jt(5),jt(6),jt(7));
        writeline(u,cmd,"127.0.0.1",7755);
    end

    function ptp(jts, vel)
        if nargin < 2
            vel = .4;
        end
        start = queryJoints;
        wayPoints = [start',jts'];
        Freq = 100;
        r = rateControl(Freq);
        T = max(abs(jts - start) / vel);
        numSamples = round(T * Freq) + 1;
        jt = trapveltraj(wayPoints,numSamples);
        for i = 1 : numSamples
            setJoints(jt(:,i));
            waitfor(r);
        end
    end
    function joints = queryJoints
        % ";" 表示查询关节角
        writeline(u,"robot;","127.0.0.1",7755);
        s = readline(u);
        joints = sscanf(s,'%f;%f;%f;%f;%f;%f;%f;')';
        joints = mod(joints + pi, 2*pi) - pi;
    end
    function T=A(r)  % r:n*1 vector
        si=size(r,1);
        modulo=norm(r);
        if modulo==0
            T=eye(si);
        else
        T=eye(si)-(1-cos(modulo))/(modulo^2)*SkewMatrix(r)+(modulo-sin(modulo))/(modulo^3)*SkewMatrix(r)*SkewMatrix(r);
        end
    end
    function vec=sym2vec(sym)
        vec=zeros(3,1);
        vec(1)=sym(3,2);
        vec(2)=sym(1,3);
        vec(3)=sym(2,1);
        
    end

function yd = impendance(t,y,M,D,K,F)



yd = zeros(6,1);
yd(1:3)=y(4:6);  
yd(4:6)=M^(-1)* (F-K*y(1:3)-D*y(4:6));
end

end