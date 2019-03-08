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
selectController = 4;
param.selectController = selectController;

%Measure matrices
Mx = zeros(2,8);
Mx(1,1) = 1; %measure of x
Mx(2,3) = 1; %measure of y
param.Mx = Mx;


Md = eye(2);
param.Md = Md;
param.dStart = [0;0];

%Tolerance
param.eps_r = eps_r;
param.eps_t = eps_t;
param.toleranceInput = 0.00002;
param.closeToTarget = 0.005;
param.angCond = 0;


param.K = [1, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 1, 0, 0, 0, 0, 0];
N =20;
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

if(selectController == 1)
    param.LTR_obsv = L_LTR_tilde;%It is not used
    
end
if(selectController ==2)
    
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


m23 = (y2 - y3)/(x2 - x3);

%-----------------using middle point

% param.xTarSigned = (x5 + 1.1*x2)*0.5;
% param.yTarSigned = (y5 + 1.1*y2)*0.5;

%-----------------using intersection 
if(x2 == x5)
    pointTarSigned = [x2 ; x2*m23-m23*xTar+yTar];
    disp('vertical');
    
else
    m25 = (y2-y5)/(x2-x5);
    ATarSigned = [[-m25 1];[-m23 1]];
    bTarSigned = [-m25*x5 + y5 ; -m23*xTar + yTar];

    pointTarSigned = ATarSigned\bTarSigned;

end



param.xTarSigned = pointTarSigned(1);
param.yTarSigned = pointTarSigned(2);
param.epsilonTarget = 0.007;



if(x1 ~= x2 && x1~=x4)

    m16 = (y1-y6)/(x1-x6);
    m21 = (y2-y1)/(x2-x1);
    
    %first rect
    c1Lower = -x1*m16+y1;
    c1Upper = -x2*m16+y2;

    c2Lower = -x3*m21+y3;
    c2Upper = -x1*m21+y1;
    
    %second rect
    c1Lower2 = -x4*m16+y4;
    c1Upper2 = -x3*m16+y3;
    
    c2Lower2 = -x4*m21+y4;
    c2Upper2 = -x2*m21+y2;



    D = zeros(6,8);
    
    D(1,1) = -m16;
    D(1,3) = 1;

    D(2,1) = -m21;
    D(2,3) = 1;
    
    D(3,1) = -m16;
    D(3,3) = 1;
    D(3,5) = -len*m16;
    D(3,7) = len;
    
    D(4,1) = -m21;
    D(4,3) = 1;
    D(4,5) = -len*m21;
    D(4,7) = len;
  

    
    

end
D(5,5) = 1;
D(6,7) = 1;

param.D = D;
D2 = D;

angleDegreeConst = deg2rad(1);
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
Q(1,1) = 1;
Q(2,2) = 1;
Q(3,3) = 1;
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
r = zeros(10,1);


eps_t = 2*param.eps_t/sqrt(2);
persistent changed

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Just change the target
if(param.selectController == 1 || param.selectController == 2 || param.selectController == 3)
    %nothing
    if(isempty(changed))
        changed = 0;
    end
    
    dist = sqrt((xHat(1)-param.xTarSigned)^2 + (xHat(3)-param.yTarSigned)^2);
    fprintf('Distance :%f \n', dist);
    if(changed == 0)
        if(dist < param.epsilonTarget)
            r(1) = param.xTar;
            r(3) = param.yTar;
            changed = 1;
        else
            r(1) = param.xTarSigned;
            r(3) = param.yTarSigned;
        end
    else
        r(1) = param.xTar;
        r(3) = param.yTar;
        
    end
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Estimation with quad prog
if(param.selectController == 4)
    if(isempty(changed))
        changed = 0;
    end
    
    dist = sqrt((xHat(1)-param.xTarSigned)^2 + (xHat(3)-param.yTarSigned)^2);
    
    %calculate the xTar
    if(changed == 0)
        if(dist < param.epsilonTarget)
            xTar = param.xTar;
            yTar = param.yTar;
            changed = 1;
        else
            xTar = param.xTarSigned;
            yTar = param.yTarSigned;
        end
    else
        xTar = param.xTar;
        yTar = param.yTar;
        
    end
    
    %Estimation with disturbance rejection
    
    Aeq = [[eye(8) - param.A , -param.B];[param.Mx , zeros(2,2)]];
    beq = [param.Bd * dHat ; [xTar ; yTar] - param.Md*dHat ];
    
    Aineq = [-eye(10) ; eye(10)];
    low = [-xTar+eps_t 0 -yTar+eps_t 0 param.angCond 0 param.angCond 0 -param.ul']';
    high = [xTar+eps_t param.eps_r yTar+eps_t param.eps_r param.angCond param.eps_r param.angCond param.eps_r param.uh']';
    bineq = [low+[param.Cd*dHat ; 0 ; 0] ; high-[param.Cd*dHat ; 0 ; 0] ];
    
    [r,~,flag] = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],[]);
    
    
%     if(flag == -2)
%         [r,~,flag] = quadprog(H,f,Aineq,bineq,[],[],[],[],[]);
%         if(flag == -2)
%             Aineq = [zeros(10,8) , [zeros(8,2) ; -eye(2)] ; [zeros(10,8) , [zeros(8,2) ; eye(2)]]];
%             bineq = [zeros(8,1) ; -param.ul ; zeros(8,1) ; param.uh];
%             [r,~,~] = quadprog(H,f,Aineq,bineq,[],[],[],[],[],[]);
%         end
%         
%     end
    

    formatSpec = 'Target xT =%f | yT =%f | uxT =%f | uyT =%f \n';
    fprintf(formatSpec,r(1),r(3),r(9),r(10));
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


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(param.selectController == 1)
    x_hat(1:8) = y;
    x_hat(9:10) = zeros(2,1);
    
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if(param.selectController == 2)
    if(isempty(check))
        %first run , we use the xHatPrev
        lhs = y-param.C*param.xStart;
        
        xExtNext = param.A*param.xStart + param.B*u  + gainObsv*lhs ;
        xHatNext = xExtNext(1:8);
        
        check = 1;
    else
        lhs = y-param.C*xHatPrevPersistent;
        xExtNext = param.A*xHatPrevPersistent + param.B*u + gainObsv*lhs ;
        xHatNext = xExtNext(1:8);
    end
    
    x_hat(1:8) = xHatNext;
    x_hat(9:10) = zeros(2,1);
    
       
    xHatPrevPersistent = xHatNext;
    
end


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(param.selectController == 3 || param.selectController == 4)
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
    

end


    formatSpec = 'Estimator xHat =%f | yHat =%f | dtx =%f | dty =%f \n';
    fprintf(formatSpec,x_hat(1),x_hat(3),x_hat(9),x_hat(10));

end % End of myStateEstimator




%% Modify the following function for your controller
function u = myMPController(r, x_hat, param)

opt = mpcqpsolverOptions;
opt.MaxIter = 200;
opt.IntegrityChecks = false;%% for code generation
opt.FeasibilityTol = 1e-3;
opt.DataType = 'double';

dHat = x_hat(9:10);
currentX = x_hat(1:8) + param.Cd *dHat;

uCurrentTarget = r(9:10);
xCurrentTarget = r(1:8);

ul = param.ul + uCurrentTarget;
uh = param.uh + uCurrentTarget;

persistent checkChangeRect
persistent checkCloseTar 

if(isempty(checkCloseTar))
   checkCloseTar = 0; 
end

if(isempty(checkChangeRect))
   checkChangeRect = 0; 
end

dist = sqrt((x_hat(1)-param.xTarSigned)^2 + (x_hat(3)-param.yTarSigned)^2);
xStart = zeros(8,1);
if(checkChangeRect == 0)
    if(dist < param.epsilonTarget)
        %Second rect constr
      
        %Recalculate the constraits
        [Dt2,Et2,bt2]=genStageConstraints(param.A,param.B,param.D2,param.cl2,param.ch2,ul,uh);
        [DD2,EE2,newbb]=genTrajectoryConstraints(Dt2,Et2,bt2,param.N);
        [F,J,L]=genConstraintMatrices(DD2,EE2,param.Gamma,param.Phi,param.N);
%         
%         F = param.F2;
%         J = param.J2;
%         L = param.L2;
%         EE = param.EE2;
%         bb = param.bb2;
%         newbb = bb - EE*kron(ones(param.N,1),uCurrentTarget);

        xStart(1) = param.xTarSigned;
        xStart(3) = param.yTarSigned;
        checkChangeRect = 1;
    else
        %First rect constr
        
        %Recalculate the constraits 
        [Dt,Et,bt]=genStageConstraints(param.A,param.B,param.D,param.cl1,param.ch1,ul,uh);
        [DD,EE,newbb]=genTrajectoryConstraints(Dt,Et,bt,param.N);
        [F,J,L]=genConstraintMatrices(DD,EE,param.Gamma,param.Phi,param.N);
        
%         F = param.F;
%         J = param.J;
%         L = param.L;
%         bb = param.bb;
%         EE = param.EE;
%         newbb = bb - EE*kron(ones(param.N,1),uCurrentTarget);

        xStart = param.xStart;
    end
else
        %Second rect constr
        [Dt2,Et2,bt2]=genStageConstraints(param.A,param.B,param.D2,param.cl2,param.ch2,ul,uh);
        [DD2,EE2,newbb]=genTrajectoryConstraints(Dt2,Et2,bt2,param.N);
        [F,J,L]=genConstraintMatrices(DD2,EE2,param.Gamma,param.Phi,param.N);

%         F = param.F2;
%         J = param.J2;
%         L = param.L2;
%         EE = param.EE2;
%         bb = param.bb2;
%         newbb = bb - kron(ones(param.N,1),uCurrentTarget);

    xStart(1) = param.xTarSigned;
    xStart(3) = param.yTarSigned;

end


f =  (param.G * (currentX-xCurrentTarget)); %linear term must be a column vector
RHS = newbb+ L*xCurrentTarget+  J*xStart; %RHS of inequality
iA = false(size(newbb));
[U,~,~]=mpcqpsolver(param.H,f,-F,-RHS,[],zeros(0,1),iA,opt);

u =U(1:param.m,:);



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    %EXTRA
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                %saturation of the input
% if(abs(u(1)) < param.toleranceInput)
%    u(1) = 0; 
% end
% if(abs(u(2)) < param.toleranceInput)
%    u(2) = 0; 
% end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                %Stop do everything close to the Tar
% dist = sqrt((currentX(1)-param.xTar)^2 + (currentX(3)-param.yTar)^2);
% if(dist < param.closeToTarget || checkCloseTar == 1)
%    u = zeros(2,1);%stop do everything 
%    checkCloseTar = 1;
% end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





formatSpec = 'MPC Controller ux=%f | uy=%f\n';
fprintf(formatSpec,u(1),u(2));
fprintf('------------------------------\n');
end % End of myMPController
