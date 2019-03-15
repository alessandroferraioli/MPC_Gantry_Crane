function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
load CraneParameters;
Ts = 1/20;

param.Ts = Ts;
[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);

param.A = A;
param.B = B;
param.C = C;
%parameters 
len = 0.47;
param.m = size(B,2);

N =30;
param.N = N;

param.K = [1, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 1, 0, 0, 0, 0, 0];

xTar = targetPoint(1);
yTar = targetPoint(2);
param.xTar = targetPoint(1);
param.yTar = targetPoint(2);

%% ++++++++++++++++++SELECT DIFFERENT  CONTROLLER+++++++++++++++++++++
%1 nothing
%2 state estimator
%3 state estimator + disturbance estimator
%4 state estimator + disturbance estimator + target calculator
%5 state estimator + target calculator
%6 LQR with just change of target
selectController = 1;
param.selectController = selectController;
%LQR BackupController
param.backupController = 0;


%% +++++++++++++++++++ MATRICES OF MEASURE ++++++++++++++++++++++++++++++
%Measure matrices
Mx = zeros(2,8);
Mx(1,1) = 1; %measure of x
Mx(2,3) = 1; %measure of y
param.Mx = Mx;


Md = eye(2);
param.Md = Md;
param.dStart = [0;0];

%% +++++++++++++++++++TOLERANCE++++++++++++++++++++++++++
param.eps_r = eps_r;
param.eps_t = eps_t;
param.toleranceInput = 0.00002;
param.epsilonTarget = 0.002;%Change target

param.closeToTarget = 0.01;
param.angCond = 0;

param.tolerance = 10^-2;




%% ++++++++++++++++++  POINTS ++++++++++++++++++++++++++++++++++
xTarget = zeros(8,1);
xTarget(1,1) = targetPoint(1);
xTarget(3,1) = targetPoint(2);
param.xTarget = xTarget;

%start point 
xStart = zeros(8,1);
xStart (1,1) = startingPoint(1);
xStart (3,1) = startingPoint(2);
param.xStart = xStart;

%% ++++++++++++++++++DISTURBANCE DYNAMICS+++++++++++++++++++++++++++++++
Cd = zeros(8,2);
Cd(1,1) = 1; %disturbance on x
Cd(3,2) = 1; %disturbance on y;
param.Cd  = Cd;

Bd = zeros(8,2);
Bd(1,1) = 1;
Bd(3,2) = 1;
param.Bd  = Bd;

Md = eye(2);
param.Md = Md;
param.dStart = [0;0];


%% +++++++++++++++++++ AUGMENTED SYSTEM ++++++++++++++++++++++++
Atilde = [[A , Bd];[zeros(2,8) , eye(2)]];
Ctilde = [C , Cd];
Btilde = [B ; zeros(2)];
addB = zeros(size(Btilde,1), size(Ctilde,2) -size(Btilde,2));


disp('rank obsv (At,Ct)');
disp(rank(obsv(Atilde,Ctilde)));

param.Atilde = Atilde;
param.Btilde = Btilde;
param.Ctilde = Ctilde;


%% +++++++++++++++++++ OBSERVATOR ++++++++++++++++++++++++++

sigma = 10^4;
weightLTR=eye(8);
Wx = sigma* (B)* (B');
Wd = sigma * Bd' * (Bd);
W = sigma * Btilde * Btilde';
weight = eye(8);

eigDes = eig(Atilde);
eigDes(1) = 0.01;
eigDes(2) = 0.02;
eigDes(9) = 0.03;
eigDes(10)= 0.04;
L1 = dlqr(A',C' ,Wx , weightLTR)';
L2 = dlqr(eye(2) , Cd' , Wd , weightLTR)';

L_LTR_tilde = dlqr(Atilde',Ctilde', eye(10) , weight)';
L_LTR= [L1 ; L2];
L_place = place(Atilde',Ctilde',eigDes)';

if(selectController == 1 || selectController == 6)
    param.LTR_obsv = L_LTR_tilde;%It is not used
end
if(selectController ==2 || selectController ==5 )
    
    param.LTR_obsv = L1;
end
if(selectController == 3 || selectController == 4)
    param.LTR_obsv = L_LTR_tilde;
end




disp('Eigen value place');
disp(eig(Atilde' - L_place*Ctilde));

disp('Eigen value LTR tilde');
disp(eig(Atilde' - L_LTR_tilde*Ctilde));

disp('Eigen value LTR ');
disp(eig(Atilde' - L_LTR*Ctilde));

%% +++++++++++++++++++input constraints++++++++++++++++++++++++++
inputConst  = 0.9 ;
ul=[-inputConst;-inputConst];
uh=[inputConst;inputConst];

param.ul = ul;
param.uh = uh;


%% ++++++++++++++++++++++++++++++STATE CONSTRAINTS ++++++++++++++++++++++
x_vec = c(:,1);
y_vec = c(:,2);
x = @(i) x_vec(i);
y = @(i) y_vec(i);

x5 = c(5,1);
y5 = c(5,2);

x2 = c(2,1);
y2 = c(2,2);

k = @(i,j) (y(j)-y(i))/(x(j)-x(i));
b = @(i,j) (x(j)*y(i)-x(i)*y(j))/(x(j)-x(i));

if (x(1) == x(2)) || (y(1) == y(2))
    'aligned'
    D1 = [1 0 0 0 0 0 0 0;
           0 0 1 0 0 0 0 0;
           1 0 0 0 r 0 0 0;
           0 0 1 0 0 0 r 0];
       
    D2 = [1 0 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0;
          1 0 0 0 r 0 0 0;
          0 0 1 0 0 0 r 0];
      
    cl1 = [min(x(6),x(2)) ; min(y(6),y(2)) ; min(x(6),x(2)) ; min(y(6),y(2))]
    ch1 = [max(x(6),x(2)) ; max(y(6),y(2)) ; max(x(6),x(2)) ; max(y(6),y(2))] 
    
    cl2 = [min(x(4),x(2)) ; min(y(4),y(2)) ; min(x(4),x(2)) ; min(y(4),y(2))]
    ch2 = [max(x(4),x(2)) ; max(y(4),y(2)) ; max(x(4),x(2)) ; max(y(4),y(2))] 
       
else    
    D1 = [-k(1,2) 0 1 0 0 0 0 0;
          -k(2,3) 0 1 0 0 0 0 0;
          -k(1,2) 0 1 0 -k(1,2)*r 0 r 0;
          -k(2,3) 0 1 0 -k(2,3)*r 0 r 0];

    cl1 = [min(b(5,6),b(1,2)) ; min(b(6,1),b(2,3)) ; min(b(5,6),b(1,2)) ; min(b(6,1),b(2,3))];
    ch1 = [max(b(5,6),b(1,2)) ; max(b(6,1),b(2,3)) ; max(b(5,6),b(1,2)) ; max(b(6,1),b(2,3))];

    D2 = [-k(1,2) 0 1 0 0 0 0 0;
          -k(2,3) 0 1 0 0 0 0 0;
          -k(1,2) 0 1 0 -k(1,2)*r 0 r 0;
          -k(2,3) 0 1 0 -k(2,3)*r 0 r 0];

    cl2 = [min(b(3,4),b(1,2)) ; min(b(4,5),b(2,3)) ; min(b(3,4),b(1,2)) ; min(b(4,5),b(2,3))];
    ch2 = [max(b(3,4),b(1,2)) ; max(b(4,5),b(2,3)) ; max(b(3,4),b(1,2)) ; max(b(4,5),b(2,3))];
end


D_ang = [0 0 0 0 1 0 0 0;  
         0 0 0 0 0 0 1 0];
     
angleDegreeConst = deg2rad(1);
cl_ang = [-angleDegreeConst ; -angleDegreeConst];
ch_ang = [angleDegreeConst ; angleDegreeConst];
% Add the angle constr
D = [D1 ; D_ang];
D2 = [D2 ; D_ang];
cl1 = [cl1 ; cl_ang];
cl2 = [cl2 ; cl_ang];
ch1 = [ch1 ; ch_ang];
ch2 = [ch2 ; ch_ang];

param.cl1 = cl1;
param.ch1 = ch1;
param.cl2 = cl2;
param.ch2 = ch2;
param.D = D;
param.D2 = D2;


%%  CALCULATING MATRICES OF CONSTRANTS AND PENALITY

% Declare penalty matrices and tune them here:
Q=eye(8);
% Q(1,1) = 3;
% Q(2,2) = 1;
% Q(3,3) = 3;
% Q(2,2) = 1;

Q(1,1) = 10;
Q(3,3) = 10;
Q(5,5) = 5;
Q(7,7) = 5;

weightInput = 0.5;
R=eye(2)*weightInput; 
P=Q; % terminal weight


% First rect
[Dt,Et,bt]=genStageConstraints(A,B,D,cl1,ch1,ul,uh);
[DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);
[Gamma,Phi] = genPrediction(A,B,N);
[F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);       
H = chol(H,'lower');
H=(H'\eye(size(H)))';

param.Gamma = Gamma;
param.Phi = Phi;


% second rect
[Dt2,Et2,bt2]=genStageConstraints(A,B,D2,cl2,ch2,ul,uh);
[DD2,EE2,bb2]=genTrajectoryConstraints(Dt2,Et2,bt2,N);
[F2,J2,L2]=genConstraintMatrices(DD2,EE2,Gamma,Phi,N);

%Saving in parameters
param.H = H;
param.G = G;
param.F = F;
param.J = J;
param.L = L;
param.Gamma = Gamma;
param.Phi = Phi;
param.EE = EE;
param.bb = bb;



param.F2 = F2;
param.J2 = J2;
param.L2 = L2;
param.bb2 = bb2;
param.EE2 = EE2;






%% CHECK RECTANGLE

param.xTarSigned = (x5 + x2)*0.5;
param.yTarSigned = (y5 + y2)*0.5;



start = [startingPoint(1) 0 startingPoint(2) 0 0 0 0 0]';
target = [targetPoint(1) 0 targetPoint(2) 0 0 0 0 0]';

Dt1 = Dt;
bt1 = bt;
if ( (Dt1*start <= bt1) & (Dt1*target <= bt1) )
    '1 -> 1'
    param.whichRectTarget = 1;
    
    param.xTarSigned = xTar;
    param.yTarSigned = yTar;

elseif (Dt2*start <= bt2) & (Dt2*target <= bt2)
    '2 -> 2'
    
      param.whichRectTarget = 1;
    
    param.xTarSigned = xTar;
    param.yTarSigned = yTar;
elseif (Dt1*start <= bt1) & (Dt2*target <= bt2)
    '1 -> 2'
     param.whichRectTarget = 2;


elseif (Dt2*start <= bt2) & (Dt1*target <= bt1) 
    '2 -> 1'
    param.whichRectTarget = 1;    
    param.xTarSigned = xTar;
    param.yTarSigned = yTar;

end





%% BACKUP : Controller LQR close to the solution
R_backup = 0.02*eye(2);
Q_backup = C' * C;

param.K_LQR = dlqr(A,B,Q_backup,R_backup);




end % End of mySetup

%% ##################################################################################################




