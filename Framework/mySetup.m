function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
%generate matrix of dynamics
load CraneParameters;
Ts = 1/20;

param.Ts = Ts;
[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);

param.A = A;
param.B = B;
param.C = C;

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

%Measure matrices
Mx = zeros(2,8);
Mx(1,1) = 1; %measure of x
Mx(2,3) = 1; %measure of y
param.Mx = Mx;


Md = eye(2);
param.Md = Md;
param.dStart = [0;0];

%TOLERANCES
param.eps_r = eps_r;
param.eps_t = eps_t;
param.toleranceInput = 0.00002;
param.epsilonTarget = 0.002;%Change target

param.closeToTarget = 0.01;
param.angCond = 0;


param.K = [1, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 1, 0, 0, 0, 0, 0];
N =30;
param.N = N;
xTar = targetPoint(1);
yTar = targetPoint(2);
param.xTar = targetPoint(1);
param.yTar = targetPoint(2);

%parameters 
len = 0.47;
param.m = size(B,2);

%target point
xTarget = zeros(8,1);
xTarget(1,1) = targetPoint(1);
xTarget(3,1) = targetPoint(2);
param.xTarget = xTarget;

%start point 
xStart = zeros(8,1);
xStart (1,1) = startingPoint(1);
xStart (3,1) = startingPoint(2);
param.xStart = xStart;



%Disturbance
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


%tilde system
Atilde = [[A , Bd];[zeros(2,8) , eye(2)]];
Ctilde = [C , Cd];
Btilde = [B ; zeros(2)];
addB = zeros(size(Btilde,1), size(Ctilde,2) -size(Btilde,2));


disp('rank obsv (At,Ct)');
disp(rank(obsv(Atilde,Ctilde)));

param.Atilde = Atilde;
param.Btilde = Btilde;
param.Ctilde = Ctilde;


%observator

param.tolerance = 10^-2;

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

%input constraints
inputConst  = 0.8 ;
ul=[-inputConst;-inputConst];
uh=[inputConst;inputConst];

param.ul = ul;
param.uh = uh;


%state constraints
x1 = c(1,1);
y1 = c(1,2);

x2 = c(2,1);
y2 = c(2,2);

x3 = c(3,1);
y3 = c(3,2);

x4 = c(4,1);
y4 = c(4,2);

x5 = c(5,1);
y5 = c(5,2);

x6 = c(6,1);
y6 = c(6,2);

caseConst = -1;

%1-2
if(abs(x2-x5) < 10^(-10))
    if(y2>y5)
        caseConst = 1;
        disp('UP');
    else
        caseConst = 2;
        disp('DOWN');
    end
    
end


%3-4
if(abs(y2-y5) < 10^(-10))
    if(x2>x5)
        caseConst = 3;
        disp('RIGHT');
    else
        caseConst = 4;
        disp('LEFT');
    end
    
end

%5-6
if(abs(y1-y2) < 10^(-10))
    if(x2>x1)
        caseConst = 5;
        disp('UP RIGHT');
    else
        caseConst = 6;
        disp('BOTTOM LEFT');
    end
    
end
%7-8
if(abs(x1-x2) < 10^(-10))
    if(y1>y2)
        caseConst = 7;
        disp('BOTTOM RIGHT');
    else
        caseConst = 8;
        disp('UP LEFT');
    end
    
end


shrinkFactor = 0.00;%shrink factor of the constraints

if(caseConst == 1)
    
  m1 = (y2-y1)/(x2-x1);
  m2 = (y1-y6)/(x1-x6);
  
  %first rectangle
  c1Lower = (-x1*m2+y1)*(1 + shrinkFactor);
  c1Upper = (-x2*m2+y2)*(1 - shrinkFactor);
  
  c2Lower = (-x6*m1+y6)*(1 + shrinkFactor);
  c2Upper = (-x1*m1+y1)*(1 - shrinkFactor);
  
  %second rectangle
  c1Lower2 = (-x4*m2+y4)*(1 + shrinkFactor);
  c1Upper2 = (-x3*m2+y3)*(1 - shrinkFactor);
  
  c2Lower2 = (-x4*m1+y4)*(1 + shrinkFactor);
  c2Upper2 = (-x2*m1+y2)*(1 - shrinkFactor);

    
    
end


if(caseConst == 2)
    
  m1 = (y1-y2)/(x2-x1);
  m2 = (y6-y1)/(x6-x1);
  
  %first rectangle
  c1Lower = (-x2*m2+y2)*(1 + shrinkFactor);
  c1Upper = (-x1*m2+y1)*(1 - shrinkFactor);
  
  c2Lower = (-x2*m1+y2)*(1 + shrinkFactor);
  c2Upper = (-x6*m1+y6)*(1 - shrinkFactor);
  
  %second rectangle
  c1Lower2 = (-x3*m2+y3)*(1 + shrinkFactor);
  c1Upper2 = (-x4*m2+y4)*(1 - shrinkFactor);
  
  c2Lower2 = (-x2*m1+y2)*(1 + shrinkFactor);
  c2Upper2 = (-x4*m1+y4)*(1 - shrinkFactor);
    
end


if(caseConst == 3)
    
  m1 = (y1-y6)/(x1-x6);
  m2 = (y2-y1)/(x2-x1);
  
  %first rectangle
  c1Lower = (-x6*m2+y6)*(1 + shrinkFactor);
  c1Upper = (-x1*m2+y1)*(1 - shrinkFactor);
  
  c2Lower = (-x2*m1+y2)*(1 + shrinkFactor);
  c2Upper = (-x1*m1+y1)*(1 - shrinkFactor);
  
  %second rectangle
  c1Lower2 = (-x3*m2+y3)*(1 + shrinkFactor);
  c1Upper2 = (-x2*m2+y2)*(1 - shrinkFactor);
  
  c2Lower2 = (-x3*m1+y3)*(1 + shrinkFactor);
  c2Upper2 = (-x5*m1+y5)*(1 - shrinkFactor);
    
end

if(caseConst == 4)
    
  m1 = (y1-y6)/(x1-x6);
  m2 = (y2-y1)/(x2-x1);
  
  %first rectangle
  c1Lower = (-x2*m2+y2)*(1 + shrinkFactor);
  c1Upper = (-x5*m2+y5)*(1 - shrinkFactor);
  
  c2Lower = (-x6*m1+y6)*(1 + shrinkFactor);
  c2Upper = (-x2*m1+y2)*(1 - shrinkFactor);
  
  %second rectangle
  c1Lower2 = (-x2*m2+y2)*(1 + shrinkFactor);
  c1Upper2 = (-x4*m2+y4)*(1 - shrinkFactor);
  
  c2Lower2 = (-x5*m1+y5)*(1 + shrinkFactor);
  c2Upper2 = (-x3*m1+y3)*(1 - shrinkFactor);
    
end

if(caseConst == 5)
   
 %first rectangle
    c1Lower = x6*(1 + shrinkFactor);
    c1Upper = x3*(1 - shrinkFactor);
    
    c2Lower = y5*(1 + shrinkFactor);
    c2Upper = y2*(1 - shrinkFactor);
    
    %second rectangle
    c1Lower2 = x4*(1 + shrinkFactor);
    c1Upper2 = x3*(1 - shrinkFactor);
    
    c2Lower2 = y3*(1 + shrinkFactor);
    c2Upper2 = y2*(1 - shrinkFactor);
    
    D = zeros(6,8);
    D(1,1) = 1;
    D(2,3) = 1;
    D(3,1) = 1;
    D(3,5) = r;
    D(4,3) = 1;
    D(4,7) = r;
    
    
end


if(caseConst == 6)
    
 %first rectangle
    c1Lower = x2*(1 + shrinkFactor);
    c1Upper = x1*(1 - shrinkFactor);
    
    c2Lower = y1*(1 + shrinkFactor);
    c2Upper = y6*(1 - shrinkFactor);
    
    %second rectangle
    c1Lower2 = x2*(1 + shrinkFactor);
    c1Upper2 = x5*(1 - shrinkFactor);
    
    c2Lower2 = y2*(1 + shrinkFactor);
    c2Upper2 = y3*(1 - shrinkFactor);
    
    D = zeros(6,8);
    D(1,1) = 1;
    D(2,3) = 1;
    D(3,1) = 1;
    D(3,5) = r;
    D(4,3) = 1;
    D(4,7) = r;
  
end


if(caseConst == 7)
    
    %first rectangle
    c1Lower = x6*(1 + shrinkFactor);
    c1Upper = x1*(1 - shrinkFactor);
    
    c2Lower = y2*(1 + shrinkFactor);
    c2Upper = y6*(1 - shrinkFactor);
    
    %second rectangle
    c1Lower2 = x3*(1 + shrinkFactor);
    c1Upper2 = x2*(1 - shrinkFactor);
    
    c2Lower2 = y3*(1 + shrinkFactor);
    c2Upper2 = y4*(1 - shrinkFactor);
    
    D = zeros(6,8);
    D(1,1) = 1;
    D(2,3) = 1;
    D(3,1) = 1;
    D(3,5) = r;
    D(4,3) = 1;
    D(4,7) = r;
 
end


if(caseConst == 8)
    
    %first rectangle
    c1Lower = x1*(1 + shrinkFactor);
    c1Upper = x6*(1 - shrinkFactor);
    
    c2Lower = y1*(1 + shrinkFactor);
    c2Upper = y2*(1 - shrinkFactor);
    
    %second rectangle
    c1Lower2 = x2*(1 + shrinkFactor);
    c1Upper2 = x3*(1 - shrinkFactor);
    
    c2Lower2 = y4*(1 + shrinkFactor);
    c2Upper2 = y3*(1 - shrinkFactor);
    
    D = zeros(6,8);
    D(1,1) = 1;
    D(2,3) = 1;
    D(3,1) = 1;
    D(3,5) = r;
    D(4,3) = 1;
    D(4,7) = r;

end




%CHEK WHERE IS THE TARGET RESPECT TO THE RECTANGLE 


if(caseConst == 1)

    respect_rectangles = (xTar - x5)*(y5 - y6) - (yTar -y5)*(x5-x6);
    
    if(respect_rectangles >0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    
end


if(caseConst == 2)
respect_rectangles = (xTar - x5)*(y5 - y6) - (yTar -y6)*(x5-x6);
    
    if(respect_rectangles > 0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    
end


if(caseConst == 3)
    %NB we have - cause this shape is obtained by a rotation of pi/2
    respect_rectangles = -((xTar - x5)*(y5 - y6) - (yTar -y5)*(x5-x6));
    disp(respect_rectangles);
    if(respect_rectangles < 0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
  
end

if(caseConst == 4)
    respect_rectangles = -((xTar - x5)*(y5 - y6) - (yTar -y5)*(x5-x6));
    
    if(respect_rectangles <0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    

    
end

if(caseConst == 5)
    
    if(yTar < y5 && xTar>x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    
    
end


if(caseConst == 6)
    if(yTar > y5 && xTar<x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end

  
end


if(caseConst == 7)
    
    if(yTar < y5 && xTar<x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
 
end


if(caseConst == 8)
    
    if(yTar > y5 && xTar>x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end

end





% %-----------------using middle point
% 
% param.xTarSigned = (x5 + x2)*0.5;
% param.yTarSigned = (y5 + y2)*0.5;

%-----------------using intersection 

%m23 = (y2 - y3)/(x2 - x3);

% if(x2 == x5)
%     pointTarSigned = [x2 ; x2*m23-m23*xTar+yTar];
%     disp('vertical');
%     
% else
%     m25 = (y2-y5)/(x2-x5);
%     ATarSigned = [[-m25 1];[-m23 1]];
%     bTarSigned = [-m25*x5 + y5 ; -m23*xTar + yTar];
% 
%     pointTarSigned = ATarSigned\bTarSigned;
% 
% end
% param.xTarSigned = pointTarSigned(1);
% param.yTarSigned = pointTarSigned(2);


if(caseConst == 1 || caseConst  == 2|| caseConst  == 3|| caseConst  == 4)

    D = zeros(6,8);
    
    D(1,1) = -m2;
    D(1,3) = 1;
    
    D(2,1) = -m1;
    D(2,3) = 1;
    
    D(3,1) = -m2;
    D(3,3) = 1;
    D(3,5) = -len*m2;
    D(3,7) = len;
    
    D(4,1) = -m1;
    D(4,3) = 1;
    D(4,5) = -len*m1;
    D(4,7) = len;


end

    
    
D(5,5) = 1;
D(6,7) = 1;
D2 = D;

angleDegreeConst = deg2rad(1);
angleConstraint = [-angleDegreeConst , angleDegreeConst];


cl1 = [c1Lower;c2Lower;c1Lower;c2Lower;angleConstraint(1);angleConstraint(1)];
ch1 = [c1Upper;c2Upper;c1Upper;c2Upper;angleConstraint(2);angleConstraint(2)];
param.cl1 = cl1;
param.ch1 = ch1;
param.D = D;

%1-2 chart
%3-4 mass
%5-6 angle
cl2 = [c1Lower2;c2Lower2;c1Lower2;c2Lower2;angleConstraint(1);angleConstraint(1)];
ch2 = [c1Upper2;c2Upper2;c1Upper2;c2Upper2;angleConstraint(2);angleConstraint(2)];
param.cl2 = cl2;
param.ch2 = ch2;
param.D2 = D2;




%+++++++++++++++++++++++++++++++++++++++++++matrix 1
% Declare penalty matrices and tune them here:
Q=C'*C;
Q(1,1) = 3;
Q(2,2) = 1;
Q(3,3) = 3;
Q(2,2) = 1;

weightInput = 0.01;
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




%BACKUP : Controller LQR close to the solution
R_backup = 0.02*eye(2);
Q_backup = C' * C;

param.K_LQR = dlqr(A,B,Q_backup,R_backup);




end % End of mySetup

%% ##################################################################################################




