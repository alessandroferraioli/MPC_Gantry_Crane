function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
%generate matrix of dynamics
load CraneParameters;
Ts = 1/20;

param.Ts = Ts;
[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);

param.A = A;
param.B = B;
param.C = C;

%Measure matrices
Mx = zeros(2,8);
Mx(1,1) = 1; %measure of x
Mx(2,3) = 1; %measure of y
param.Mx = Mx;


Md = eye(2);
param.Md = Md;
param.dStart = [0;0];

param.eps_r = eps_r;
param.eps_t = eps_t;
param.angCond = 0;


param.K = [1, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 1, 0, 0, 0, 0, 0];
N =20;
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

%half Point
changePoint = 4;
param.HalfPoint = [(xTarget(1,1) - xStart(1,1))/changePoint , (xTarget(3,1) - xStart(3,1))/changePoint];

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


%observator

param.tolerance = 10^-2;

sigma = 10^4;
weightLTR=eye(8);
Wx = sigma* (B)* (B');
Wd = sigma * Bd' * (Bd);

L1 = dlqr(A',C' ,Wx , weightLTR)';
L2 = dlqr(eye(2) , Cd' , Wd , weightLTR)';

param.LTR_obsv = [L1 ; L2];


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

%calculation of some point of interest
m56 = (y5 - y6)/(x5 - x6);
m23 = (y2 - y3)/(x2 - x3);

m21 = (y2-y1)/(x2-x1);
m54 = (y5-y4)/(x5-x4);




Apoint = [[-m56 ,1];[-m23 , 1]];
bpoint = [-x6*m56+y6, -x3*m23 + y3]';

pointInt = Apoint\bpoint;
x7 = pointInt(1);
y7 = pointInt(2);


Apoint = [[-m21 ,1];[-m54 , 1]];
bpoint = [-x1*m21+y1, -x4*m54 + y4]';

pointInt = Apoint\bpoint;
x8 = pointInt(1);
y8 = pointInt(2);


param.xTarSigned = (x5 + 1.1*x2)*0.5;
param.yTarSigned = (y5 + 1.1*y2)*0.5;
param.epsilonTarget = 0.03;


if(x1 == x2)%vertical
    c1Lower = x1;
    c1Upper = x4;
    c2Lower = y1;
    c2Upper = y2;


    D = zeros(6,8);
    
    D(1,1) = 1;
    
    D(2,3) = 1;
    
    D(3,1) = 1;
    D(3,5) = len;
    
    D(4,3) = 1;
    D(4,7) = len;

end

if(x1 == x4)%horizontal
    c1Lower = x1;
    c1Upper = x2;
    c2Lower = y4;
    c2Upper = y1;
    
    D = zeros(6,8);
    
    D(1,1) = 1;
    
    D(2,3) = 1;
    
    D(3,1) = 1;
    D(3,5) = len;
    
    D(4,3) = 1;
    D(4,7) = len;

end

if(x1 ~= x2 && x1~=x4)
    %normal
    
    m1 = (y6-y1)/(x6-x1);
    m2 = (y1-y2)/(x1-x2);

    c1Lower = -x1*m1+y1;
    c2Lower = -x7*m2+y7;

    c1Upper = -x2*m1+y2;
    c2Upper = -x2*m2+y2;


    %matrix of state constrait

    D = zeros(6,8);
    
    D(1,1) = -m1;
    D(1,3) = 1;

    D(2,1) = -m2;
    D(2,3) = 1;
    
    D(3,1) = -m1;
    D(3,3) = 1;
    D(3,5) = -len*m1;
    D(3,7) = len;
    
    D(4,1) = -m2;
    D(4,3) = 1;
    D(4,5) = -len*m2;
    D(4,7) = len;
    
    
    
     %second rect

    
    m3 = (y4-y8)/(x4-x8);
    m4 = (y1-y2)/(x1-x2);

    c1Lower2 = -x8*m3+y8;
    c2Lower2 = -x3*m4+y3;

    c1Upper2 = -x2*m3+y2;
    c2Upper2 = -x2*m4+y2;


    %matrix of state constrait

    D2 = zeros(6,8);
    
    D2(1,1) = -m3;
    D2(1,3) = 1;

    D2(2,1) = -m4;
    D2(2,3) = 1;
    
    D2(3,1) = -m3;
    D2(3,3) = 1;
    D2(3,5) = -len*m3;
    D2(3,7) = len;
    
    D2(4,1) = -m4;
    D2(4,3) = 1;
    D2(4,5) = -len*m4;
    D2(4,7) = len;
    



end


%angle constrait
D(5,5) = 1;
D(6,7) = 1;

angleDegreeConst = deg2rad(2);
angleConstraint = [-angleDegreeConst , angleDegreeConst];
cl1 = [c1Lower;c2Lower;c1Lower;c2Lower;angleConstraint(1);angleConstraint(1)];
ch1 = [c1Upper;c2Upper;c1Upper;c2Upper;angleConstraint(2);angleConstraint(2)];


param.cl1 = cl1;
param.ch1 = ch1;

param.D1 = D;


cl2 = [c1Lower2;c2Lower2;c1Lower2;c2Lower2;angleConstraint(1);angleConstraint(1)];
ch2 = [c1Upper2;c2Upper2;c1Upper2;c2Upper2;angleConstraint(2);angleConstraint(2)];


param.cl2 = cl2;
param.ch2 = ch2;

param.D2 = D2;




%+++++++++++++++++++++++++++++++++++++++++++matrix 1
% Declare penalty matrices and tune them here:
Q=C'*C;
Q(1,1) = 5;
Q(3,3) = 5;


weightInput = 0.02;
R=eye(2)*weightInput; 
P=Q; % terminal weight


% Compute stage constraint matrices and vector
[Dt,Et,bt]=genStageConstraints(A,B,D,cl1,ch1,ul,uh);
[DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);
[Gamma,Phi] = genPrediction(A,B,N);
[F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);       
H = chol(H,'lower');
H=(H'\eye(size(H)))';


%second rect
[Dt2,Et2,bt2]=genStageConstraints(A,B,D2,cl2,ch2,ul,uh);
[DD2,EE2,bb2]=genTrajectoryConstraints(Dt2,Et2,bt2,N);
[Gamma,Phi] = genPrediction(A,B,N);
[F2,J2,L2]=genConstraintMatrices(DD2,EE2,Gamma,Phi,N);
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);       
H = chol(H,'lower');
H=(H'\eye(size(H)))';


param.H = H;
param.G = G;
param.F = F;
param.J = J;
param.L = L;
param.Gamma = Gamma;
param.Phi = Phi;
param.bb = bb;




param.F2 = F2;
param.J2 = J2;
param.L2 = L2;
param.bb2 = bb2;





end % End of mySetup

%% ##################################################################################################




function r = myTargetGenerator(x_hat, param)
%% Modify the following function for your target generation
dHat = x_hat(9:10);
xHat = x_hat(1:8);

H = zeros(10,10);
H(9,9) = 1; %min on ux
H(10,10) = 1; %min on uy

f = zeros(10,1);
persistent changed
if(isempty(changed))
    
   changed = 0; 
end

dist = sqrt((xHat(1)-param.xTarSigned)^2 + (xHat(3)-param.yTarSigned)^2);
fprintf('Distance :%f \n', dist);
if(changed == 0)
    if(dist < param.epsilonTarget)
        xTar = param.xTar;
        yTar = param.yTar;
        fprintf('Final Target');
        changed = 1;
    else
        xTar = param.xTarSigned;
        yTar = param.yTarSigned;
        fprintf('Middle Target Target');
        
    end
    
else
        xTar = param.xTar;
        yTar = param.yTar;
    fprintf('Final Target');
    
end



Aeq = [[eye(8) - param.A , -param.B];[param.Mx , zeros(2,2)]];
beq = [param.Bd * dHat ; [xTar ; yTar] - param.Md*dHat ];


Aineq = [-eye(10) ; eye(10)];


%x xdot y ydot theta thetaDot phi phiDot 
%x and y epsilon_t
% x and y payload epsilon_t
%velocity x/y and velocity theta/phi epsilon_r
%input epsilon_r

low = [-xTar+param.eps_t 0 -yTar+param.eps_t 0 param.angCond 0 param.angCond 0 -param.ul']';
high = [xTar+param.eps_t param.eps_r yTar+param.eps_t param.eps_r param.angCond param.eps_r param.angCond param.eps_r param.uh']';
bineq = [low+[param.Cd*dHat ; 0 ; 0] ; high-[param.Cd*dHat ; 0 ; 0] ];
[r,~,flag] = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],[]);


if(flag == -2)
    [r,~,flag] = quadprog(H,f,Aineq,bineq,[],[],[],[],[]);
    if(flag == -2)
        Aineq = [[zeros(10,8) , [zeros(8,2) ; -eye(2)] ; [zeros(10,8) , [zeros(8,2) ; eye(2)]]]];
        bineq = [zeros(8,1) ; -param.ul ; zeros(8,1) ; param.uh];
        [r,~,flag] = quadprog(H,f,Aineq,bineq,[],[],[],[],[],[]);

    end
    
end


end % End of myTargetGenerator





function x_hat = myStateEstimator(u, y, param)

%%
persistent xHatPrevPersistent;
persistent dHatPrevPersistent;
persistent check;

x_hat = zeros(16,1);
x_hat( 1:length(y),1 ) = y; %first 8 component the measure
xExtNext = zeros(10,1);


gainObsv =  param.LTR_obsv;



if(isempty(check))
    %first run , we use the xHatPrev 
    lhs = y-param.C*param.xStart - param.Cd * param.dStart;
    
    xExtNext = [param.A*param.xStart + param.B*u + param.Bd*param.dStart; param.dStart] + gainObsv*lhs ; 
    xHatNext = xExtNext(1:8);
    dHatNext = xExtNext(9:10);
    
    check = 1;
else
    lhs = y-param.C*xHatPrevPersistent - param.Cd * dHatPrevPersistent;
    xExtNext = [param.A*xHatPrevPersistent + param.B*u  + param.Bd * dHatPrevPersistent; dHatPrevPersistent] + gainObsv*lhs ; 
    xHatNext = xExtNext(1:8);
    dHatNext = xExtNext(9:10);
 end

x_hat(1:8) = xHatNext;
x_hat(9:10) = dHatNext;

%check with the tolerance
if(abs(xHatNext(1,1) - y(1,1)) > param.tolerance ||  abs(xHatNext(3,1) - y(3,1)) > param.tolerance)
    x_hat(1:8,1) = xHatNext;
else
    x_hat(1:8,1) = xHatNext;
end

xHatPrevPersistent = xHatNext;
dHatPrevPersistent = dHatNext;
end % End of myStateEstimator




%% Modify the following function for your controller
function u = myMPController(r, x_hat, param)

opt = mpcqpsolverOptions;
opt.MaxIter = 200;
opt.IntegrityChecks = false;%% for code generation
opt.FeasibilityTol = 1e-3;
opt.DataType = 'double';

currentX = x_hat(1:8);
dHat = x_hat(9:10);
uCurrentTarget = r(9:10);
xCurrentTarget = r(1:8);


formatSpec = 'xHat =%f | yHat =%f | xT =%f | yT =%f | dtx =%f | dty =%f | utx =%f | uty =%f \n';
fprintf(formatSpec,currentX(1),currentX(3),r(1),r(3) ,dHat(1) ,dHat(2) , uCurrentTarget(1),  uCurrentTarget(2))

f =  (param.G * (currentX-xCurrentTarget)); %linear term must be a column vector
RHS = param.bb+ param.L*xCurrentTarget+  param.J*param.xStart; %RHS of inequality
iA = false(size(param.bb));
[U,~,~]=mpcqpsolver(param.H,f,-param.F,-RHS,[],zeros(0,1),iA,opt);

u =U(1:param.m,:);

end % End of myMPController
