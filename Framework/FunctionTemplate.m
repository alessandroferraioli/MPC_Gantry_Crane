
function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
%%
load CraneParameters;
Ts = 1/20;


%HW OR SW 
param.isHardware = 0;
param.Ts = Ts;

% if(param.isHardware == 1)
%     Vm = Vm*0.2;
%     Tx = Tx * 0.6;
%     Ty = Ty * 0.6;
% 
% end
% 
    Vm = Vm*0.1;
    Tx = Tx * 0.1;
    Ty = Ty * 0.1;
  m = 1*m;
    


[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);

param.A = A;
param.B = B;
param.C = C;
%parameters 
len = 0.47;
param.m = size(B,2);

if(param.isHardware ==1)
    N = 40;
    param.N = N;
else
    N =30;
    param.N = N;
end
param.K = [1, 0, 0, 0, 0, 0, 0, 0;
           0, 0, 1, 0, 0, 0, 0, 0];

xTar = targetPoint(1);
yTar = targetPoint(2);
param.xTar = targetPoint(1);
param.yTar = targetPoint(2);

%% ++++++++++++++++++SELECT DIFFERENT  CONTROLLER+++++++++++++++++++++
%1 just MPC
%2 state estimator
%3 state estimator + disturbance estimator
%4 state estimator + disturbance estimator + target calculator
%5 state estimator + target calculator
%6 LQR with just change of target

selectController = 1;
param.selectController = selectController;
%LQR BackupController
param.backupController = 0;

%CLOSE TO THE TARGET
%-1 disactive 
%0 u=0
%1 LQR
%2 PD
param.controlCloseTarget =-1;




if(param.isHardware == 1);disp('HARDWARE TEST');else;disp('SOFTWARE TEST');end

if(param.backupController == 1)
    disp('USING BACKUP CONTROLLER')
    
else
    if(param.selectController == 1);disp('USING MPC');end
    if(param.selectController == 2);disp('USING MPC + STATE ESTIMATOR');end
    if(param.selectController == 3);disp('USING  MPC + STATE/DIST ESTIMATOR');end
    if(param.selectController == 4);disp('USING MPC + STATE/DIST ESTIMATOR + TARGET');end
    if(param.selectController == 5);disp('USING MPC + STATE ESTIMATOR + TARGET');end
    if(param.selectController == 6);disp('USING LQR');end
    
    if(param.controlCloseTarget == 2);disp('USING PD CLOSE TARGET');end
    if(param.controlCloseTarget == 1);disp('USING LQR CLOSE TARGET');end
    if(param.controlCloseTarget == 0);disp('DISACTIVE CONTROL CLOSE TARGET');end
    
end

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
param.toleranceInput = 0.00002;%saturation to 0 of the input

if(param.isHardware == 1)
        param.epsilonTarget = 0.02;%Change target

else
    param.epsilonTarget = 0.002;%Change target
    
end
%param.closeToTarget = param.eps_t; %change to LQR/PID/Nothing close to the target
param.closeToTarget = 0.003;
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



%% +++++++++++++++++++ AUGMENTED SYSTEM ++++++++++++++++++++++++
Atilde = [[A , Bd];[zeros(2,8) , eye(2)]];
Ctilde = [C , Cd];
Btilde = [B ; zeros(2)];


% disp('rank obsv (At,Ct)');
% disp(rank(obsv(Atilde,Ctilde)));

param.Atilde = Atilde;
param.Btilde = Btilde;
param.Ctilde = Ctilde;


%% +++++++++++++++++++ OBSERVER ++++++++++++++++++++++++++

sigma = 10^4;
weightLTR=eye(8);
Wx = sigma* (B)* (B');
Wd = sigma * Bd' * (Bd);
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

if(selectController == 1 || selectController == 6 )
    param.LTR_obsv = L_LTR_tilde;%It is not used
end
if(selectController ==2 || selectController ==5 )
    
    param.LTR_obsv = L1;
end
if(selectController == 3 || selectController == 4)
    param.LTR_obsv = L_LTR_tilde;
end



% 
% disp('Eigen value place');
% disp(eig(Atilde' - L_place*Ctilde));
% 
% disp('Eigen value LTR tilde');
% disp(eig(Atilde' - L_LTR_tilde*Ctilde));
% 
% disp('Eigen value LTR ');
% disp(eig(Atilde' - L_LTR*Ctilde));

%% +++++++++++++++++++input constraints++++++++++++++++++++++++++
inputConst  = 0.9 ;
ul=[-inputConst;-inputConst];
uh=[inputConst;inputConst];

param.ul = ul;
param.uh = uh;


%% ++++++++++++++++++++++++++++++STATE CONSTRAINTS ++++++++++++++++++++++
x_array = c(:,1);
y_array = c(:,2);
x = @(i) x_array(i);
y = @(i) y_array(i);

%used to calculate middle point
x5 = c(5,1);
y5 = c(5,2);

x2 = c(2,1);
y2 = c(2,2);

k = @(i,j) (y(j)-y(i))/(x(j)-x(i));
b = @(i,j) (x(j)*y(i)-x(i)*y(j))/(x(j)-x(i));



if (x(1) == x(2)) || (y(1) == y(2))
    'aligned'
    D = [1 0 0 0 0 0 0 0;
           0 0 1 0 0 0 0 0;
           1 0 0 0 r 0 0 0;
           0 0 1 0 0 0 r 0];
       
    D2 = [1 0 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0;
          1 0 0 0 r 0 0 0;
          0 0 1 0 0 0 r 0];
      
    cl1 = [min(x(6),x(2)); min(y(6),y(2)) ; min(x(6),x(2)) ; min(y(6),y(2))]
    ch1 = [max(x(6),x(2)) ; max(y(6),y(2)) ; max(x(6),x(2)) ; max(y(6),y(2))] 
    
    cl2 = [min(x(4),x(2))  ; min(y(4),y(2))  ; min(x(4),x(2))  ; min(y(4),y(2)) ]
    ch2 = [max(x(4),x(2)) ; max(y(4),y(2)); max(x(4),x(2)) ; max(y(4),y(2))] 
       
else    
    D = [-k(1,2) 0 1 0 0 0 0 0;
          -k(2,3) 0 1 0 0 0 0 0;
          -k(1,2) 0 1 0 -k(1,2)*r 0 r 0;
          -k(2,3) 0 1 0 -k(2,3)*r 0 r 0;
          0 0 0 0 1 0 0 0;  
         0 0 0 0 0 0 1 0];

    cl1 = [min(b(5,6),b(1,2)) ; min(b(6,1),b(2,3)) ; min(b(5,6),b(1,2)) ; min(b(6,1),b(2,3))];
    ch1 = [max(b(5,6),b(1,2)) ; max(b(6,1),b(2,3)) ; max(b(5,6),b(1,2)) ; max(b(6,1),b(2,3))];

    D2 = [-k(1,2) 0 1 0 0 0 0 0;
          -k(2,3) 0 1 0 0 0 0 0;
          -k(1,2) 0 1 0 -k(1,2)*r 0 r 0;
          -k(2,3) 0 1 0 -k(2,3)*r 0 r 0;
            0 0 0 0 1 0 0 0;  
             0 0 0 0 0 0 1 0];

    cl2 = [min(b(3,4),b(1,2)) ; min(b(4,5),b(2,3)) ; min(b(3,4),b(1,2)) ; min(b(4,5),b(2,3))];
    ch2 = [max(b(3,4),b(1,2)) ; max(b(4,5),b(2,3)) ; max(b(3,4),b(1,2)) ; max(b(4,5),b(2,3))];
end


     
angleDegreeConst = deg2rad(1);

cl_ang = [-angleDegreeConst ; -angleDegreeConst];
ch_ang = [angleDegreeConst ; angleDegreeConst];
% Add the angle constr

shrinkFactor = 0.0;
param.shrinkFactor = shrinkFactor;
cl1 = [cl1 ; cl_ang] + shrinkFactor * [1 1 1 1 0 0]';
cl2 = [cl2 ; cl_ang] + shrinkFactor * [1 1 1 1 0 0]';
ch1 = [ch1 ; ch_ang] - shrinkFactor * [1 1 1 1 0 0]';
ch2 = [ch2 ; ch_ang] - shrinkFactor * [1 1 1 1 0 0]';

param.cl1 = cl1;
param.ch1 = ch1;
param.cl2 = cl2;
param.ch2 = ch2;
param.D = D;
param.D2 = D2;


%%  CALCULATING MATRICES OF CONSTRANTS AND PENALITY

% Declare penalty matrices and tune them here:

if(param.isHardware == 1)
    Q = eye(8);
    Q(1,1) = 25;
    Q(3,3) = 25;
    Q(5,5) = 8;
    Q(7,7) = 8;    
    weightInput = 0.02;
    R=eye(2)*weightInput;
    P=Q; % terminal weight
    
else
    Q=eye(8);
    Q(1,1) = 20;
    Q(3,3) = Q(1,1);
    Q(5,5) = 15;
    Q(7,7) = Q(7,7);
    %weightInput = 0.5;
    weightInput = 1.2;
    R=eye(2)*weightInput;
    P=Q; % terminal weight
end

disp('Matrix Q : ');disp(Q);
disp('Matrix R : ');disp(R);
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



start_extended = [startingPoint(1) 0 startingPoint(2) 0 0 0 0 0]';
target_extended = [targetPoint(1) 0 targetPoint(2) 0 0 0 0 0]';


if ( (Dt*start_extended <= bt) & (Dt*target_extended <= bt) )
    disp('TARGET FIRST RECT , START FIRST RECT');
    param.whichRectTarget = 1;
    
    param.xTarSigned = xTar;
    param.yTarSigned = yTar;
    
elseif (Dt2*start_extended <= bt2) & (Dt2*target_extended <= bt2)
    disp('TARGET SECOND RECT , START SECOND RECT');
    
    param.whichRectTarget = 1;
    
    param.xTarSigned = xTar;
    param.yTarSigned = yTar;
elseif (Dt*start_extended <= bt) & (Dt2*target_extended <= bt2)
    disp('TARGET SECOND RECT , START FIRST RECT');
    
    param.whichRectTarget = 2;
    
    
elseif (Dt2*start_extended <= bt2) & (Dt*target_extended <= bt)
    disp('TARGET FIRST RECT , START SECOND RECT');
    
    param.whichRectTarget = 1;
    param.xTarSigned = xTar;
    param.yTarSigned = yTar;

end






%% BACKUP : Controller LQR close to the solution

    Q_backup=eye(8);
    Q_backup(1,1) = 50;
    Q_backup(3,3) = Q_backup(1,1);
    Q_backup(5,5) = 60;
    Q_backup(7,7) = Q_backup(7,7);
    weightInput_backup = 1.5;
    R_backup=eye(2)*weightInput_backup;


param.K_LQR = dlqr(A,B,Q_backup,R_backup);




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


eps_t = (param.eps_t);%/sqrt(2);
persistent changed
persistent timeCloseMiddle

if(isempty(timeCloseMiddle))
    timeCloseMiddle = 0;
end
dist = sqrt((xHat(1)-param.xTarSigned)^2 + (xHat(3)-param.yTarSigned)^2);
dist_final = sqrt((xHat(1)-param.xTar)^2 + (xHat(3)-param.yTar)^2);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%JUST CHANGE THE TARGET FROM THE MIDDLE TO THE END
if(param.selectController == 1 || param.selectController == 2 || param.selectController == 3 || param.selectController == 6 )
    %nothing
    if(isempty(changed))
        changed = 0;
    end
    
    fprintf('Distance Middle :%f | Distance Final :%f\n', dist,dist_final);
    if(changed == 0)
        if(dist < param.epsilonTarget )%Change target to the final
                r(1) = param.xTar;
                r(3) = param.yTar;
                changed = 1;%to be sure that i am not going to change again the target

        else
            r(1) = param.xTarSigned; %target still middle Point
            r(3) = param.yTarSigned;
            
        end
    else %if we changed at least one time --> always final target
        r(1) = param.xTar;
        r(3) = param.yTar;
        disp('FINAL TARGET');
        
    end
    
    if(param.whichRectTarget == 1)
       r(1) = param.xTar;
       r(3) = param.yTar;
    end
    
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%ESTIMATION WITH QUAD PROG
if(param.selectController == 4 || param.selectController==5)
    if(isempty(changed))
        changed = 0;
    end
    
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
    low = [xTar-eps_t 0 yTar-eps_t 0 param.angCond 0 param.angCond 0 -param.ul']';
    high = [xTar+eps_t param.eps_r yTar+eps_t param.eps_r param.angCond param.eps_r param.angCond param.eps_r param.uh']';
    bineq = [low+[param.Cd*dHat ; 0 ; 0] ; high-[param.Cd*dHat ; 0 ; 0] ];
    
    [r,~,flag] = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],[]);
    if(flag == -2)
        r(1) = xTar;
        r(3) = yTar;
        
    end
    
    if(flag == -2)
        [r,~,flag] = quadprog(H,f,Aineq,bineq,[],[],[],[],[]);
        if(flag == -2)
            Aineq = [zeros(10,8) , [zeros(8,2) ; -eye(2)] ; [zeros(10,8) , [zeros(8,2) ; eye(2)]]];
            bineq = [zeros(8,1) ; -param.ul ; zeros(8,1) ; param.uh];
            [r,~,~] = quadprog(H,f,Aineq,bineq,[],[],[],[],[],[]);
        end
        
    end
    

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


%+++++++++++++++++++++++++++++++++JUST MEASURE ++++++++++++++++++++++++++++++++++++
if(param.selectController == 1 || param.selectController == 6)
    %no estimator
    x_hat(1:8) = y;
    x_hat(9:10) = zeros(2,1);
    
end
%+++++++++++++++++++++++++++++JUST STATE ESTIMATOR++++++++++++++++++++++++++++++++++++++++

if(param.selectController == 2 || param.selectController == 5)
    if(isempty(check))
        %first run , we use the xHatStartt
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


%+++++++++++++++++++++++++++++++++++STATE AND DISTURBANCE ESTIAMTOR++++++++++++++++++++++++++++++++++
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

persistent checkChangeRect %to change target just one time 
persistent checkCloseTar  %to change controller just one time close to the targeet
persistent prevError %used in the PD 
persistent prevErrorController

if(isempty(prevErrorController))
   prevErrorController = 0; 
end

if(isempty(prevError))
   prevError = 0; 
end



if(isempty(checkCloseTar))
   checkCloseTar = 0; 
end

if(isempty(checkChangeRect))
   checkChangeRect = 0; 
end

recalculateConstr = 0;
if(param.selectController == 4 || param.selectController == 5)
    recalculateConstr = 1;        
end

dist_middle = sqrt((x_hat(1)-param.xTarSigned)^2 + (x_hat(3)-param.yTarSigned)^2);
dist_final = sqrt((x_hat(1)-param.xTar)^2 + (x_hat(3)-param.yTar)^2);

xStart = zeros(8,1);


%+++++++++++++++++++MPC++++++++++++++++++++++++++++++++++++++++++++
if(param.selectController ~= 6 )
    if(checkChangeRect == 0 )
        if(dist_middle < param.epsilonTarget && param.whichRectTarget == 2)
            %Second rect constr
            
            %Recalculate the constraits
            if(recalculateConstr == 1)
                [Dt2,Et2,bt2]=genStageConstraints(param.A,param.B,param.D2,param.cl2,param.ch2,ul,uh);
                [DD2,EE2,newbb]=genTrajectoryConstraints(Dt2,Et2,bt2,param.N);
                [F,J,L]=genConstraintMatrices(DD2,EE2,param.Gamma,param.Phi,param.N);
            end
            xStart(1) = param.xTarSigned;
            xStart(3) = param.yTarSigned;
            checkChangeRect = 1;
        else%else dist
            %First rect constr
            
            %Recalculate the constraits
            if(recalculateConstr == 1)
                [Dt,Et,bt]=genStageConstraints(param.A,param.B,param.D,param.cl1,param.ch1,ul,uh);
                [DD,EE,newbb]=genTrajectoryConstraints(Dt,Et,bt,param.N);
                [F,J,L]=genConstraintMatrices(DD,EE,param.Gamma,param.Phi,param.N);
            end
            xStart = param.xStart;
        end
    else
        %Second rect constr
        if(recalculateConstr == 1)
            [Dt2,Et2,bt2]=genStageConstraints(param.A,param.B,param.D2,param.cl2,param.ch2,ul,uh);
            [DD2,EE2,newbb]=genTrajectoryConstraints(Dt2,Et2,bt2,param.N);
            [F,J,L]=genConstraintMatrices(DD2,EE2,param.Gamma,param.Phi,param.N);
        end
        xStart(1) = param.xTarSigned;
        xStart(3) = param.yTarSigned;
        
    end
    if(param.whichRectTarget == 1)
        xStart = param.xStart;
        
    end
    
    if(recalculateConstr == 1)
        f =  (param.G * (currentX-xCurrentTarget)); %linear term must be a column vector
        RHS = newbb+ L*xCurrentTarget+  J*xStart; %RHS of inequality
        iA = false(size(newbb));
        [U,~,value]=mpcqpsolver(param.H,f,-F,-RHS,[],zeros(0,1),iA,opt);
        u =U(1:param.m,:);
    elseif(checkChangeRect == 1)
        
        f =  (param.G * (currentX-xCurrentTarget)); %linear term must be a column vector
        RHS = param.bb2+ param.L2*xCurrentTarget+  param.J2*xStart; %RHS of inequality
        iA = false(size(param.bb2));
        [U,~,value]=mpcqpsolver(param.H,f,-param.F2,-RHS,[],zeros(0,1),iA,opt);
        u =U(1:param.m,:);
    elseif(checkChangeRect == 0)
        f =  (param.G * (currentX-xCurrentTarget)); %linear term must be a column vector
        RHS = param.bb+ param.L*xCurrentTarget+  param.J*xStart; %RHS of inequality
        iA = false(size(param.bb));
        [U,~,value]=mpcqpsolver(param.H,f,-param.F,-RHS,[],zeros(0,1),iA,opt);
        u =U(1:param.m,:);
    end
end

if(param.selectController == 6)

    u = -param.K_LQR * (currentX - xCurrentTarget);
    disp('LQR - CHANGE TARGET');
    if(abs(u(1)) >= 1 )%control is just on x axis
        disp('SATURATED INPUT LQR CONTROLLER');
        u(1) = 0.8 * sign(u(1));
    end
    if(abs(u(2)) > 1 )%control is just on x axis
        disp('SATURATED INPUT LQR CONTROLLER');
        u(2) = 0.8 * sign(u(2));
    end
    
end






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
if(param.controlCloseTarget == 0) 
    if(dist_final < param.closeToTarget || checkCloseTar == 1 && checkChangeRect == 1)
        u = zeros(2,1);%stop do everything
        checkCloseTar = 1;
    end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%LQR CLOSE FINAL TARGET
if(param.controlCloseTarget == 1 )
    if(dist_final < param.closeToTarget || checkCloseTar == 1 && checkChangeRect == 1)
        disp('LQR CLOSE SOLUTION');
        u = -param.K_LQR * (currentX - param.xTarget); 
        checkCloseTar = 1;
    end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PD CONTROLLER CLOSE FINAL TARGET
Kp = 1;
Kd = 0.5;
error =  (currentX - param.xTarget);

approx_deriv_err  =(error - prevError)/param.Ts;
if(param.controlCloseTarget == 2 )
    if(dist_final < param.closeToTarget || (checkCloseTar == 1 && checkChangeRect == 1))
        %checkChangeRect--> we need to active this control iff we changed
        %target
        uPID = Kp * error + Kd * approx_deriv_err;
        u = [uPID(5);uPID(7)]; %Pd controller on error of angle
        checkCloseTar = 1;
        disp('PD CLOSE SOLUTION');
    end
end
prevError = error;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %BACKUP : Just use LQR 
if( param.backupController == 1 && checkChangeRect == 0)   
   disp('BACKUP - LQR WITHOUT CONSTRAINTS');
    u = -param.K_LQR * (currentX - param.xTarget); 
    if(abs(u(1)) >= 1 )%control is just on x axis
        disp('SATURATED INPUT BACKUP CONTROLLER');
       u(1) = 0.5 * sign(u(1)); 
    end
        if(abs(u(2)) >= 1 )%control is just on x axis
        disp('SATURATED INPUT BACKUP CONTROLLER');
       u(2) = 0.5 * sign(u(2)); 
    end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




formatSpec = 'MPC Controller ux=%f | uy=%f\n';
fprintf(formatSpec,u(1),u(2));
fprintf('------------------------------\n');
end % End of myMPController
