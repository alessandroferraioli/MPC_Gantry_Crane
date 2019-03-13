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
persistent prevError

if(isempty(prevError))
   prevError = 0; 
end

if(isempty(checkCloseTar))
   checkCloseTar = 0; 
end

if(isempty(checkChangeRect))
   checkChangeRect = 0; 
end

dist_middle = sqrt((x_hat(1)-param.xTarSigned)^2 + (x_hat(3)-param.yTarSigned)^2);
dist_final = sqrt((x_hat(1)-param.xTar)^2 + (x_hat(3)-param.yTar)^2);

xStart = zeros(8,1);
if(checkChangeRect == 0 )
    if(dist_middle < param.epsilonTarget && param.whichRectTarget == 2)
        %Second rect constr
      
        %Recalculate the constraits
        [Dt2,Et2,bt2]=genStageConstraints(param.A,param.B,param.D2,param.cl2,param.ch2,ul,uh);
        [DD2,EE2,newbb]=genTrajectoryConstraints(Dt2,Et2,bt2,param.N);
        [F,J,L]=genConstraintMatrices(DD2,EE2,param.Gamma,param.Phi,param.N);

        xStart(1) = param.xTarSigned;
        xStart(3) = param.yTarSigned;
        checkChangeRect = 1;
    else
        %First rect constr
        
        %Recalculate the constraits 
        [Dt,Et,bt]=genStageConstraints(param.A,param.B,param.D,param.cl1,param.ch1,ul,uh);
        [DD,EE,newbb]=genTrajectoryConstraints(Dt,Et,bt,param.N);
        [F,J,L]=genConstraintMatrices(DD,EE,param.Gamma,param.Phi,param.N);

        xStart = param.xStart;
    end
else
        %Second rect constr
        [Dt2,Et2,bt2]=genStageConstraints(param.A,param.B,param.D2,param.cl2,param.ch2,ul,uh);
        [DD2,EE2,newbb]=genTrajectoryConstraints(Dt2,Et2,bt2,param.N);
        [F,J,L]=genConstraintMatrices(DD2,EE2,param.Gamma,param.Phi,param.N);

        xStart(1) = param.xTarSigned;
        xStart(3) = param.yTarSigned;

end
if(param.whichRectTarget == 1)
     xStart = param.xStart;

end

if(param.selectController == 6 && param.backupController==0)
  u = -param.K_LQR * (currentX - xCurrentTarget); 
else
    f =  (param.G * (currentX-xCurrentTarget)); %linear term must be a column vector
    RHS = newbb+ L*xCurrentTarget+  J*xStart; %RHS of inequality
    iA = false(size(newbb));
    [U,~,value]=mpcqpsolver(param.H,f,-F,-RHS,[],zeros(0,1),iA,opt);
    u =U(1:param.m,:);

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
%Switch to LQR if close to the solution

if(param.controlCloseTarget == 1)
    if(dist_final < param.closeToTarget || checkCloseTar == 1 && checkChangeRect == 1)
        disp('LQR CLOSE SOLUTION');
        u = -param.K_LQR * (currentX - param.xTarget); 
        checkCloseTar = 1;
    end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%PD CONTROLLER
Kp = 0;
Kd = 1;
error =  (currentX - param.xTarget);

approx_deriv_err  =(error - prevError)/param.Ts;
if(param.controlCloseTarget == 2)
    if(dist_final < param.closeToTarget || checkCloseTar == 1 && checkChangeRect == 1)
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
if( param.backupController == 1)   
   disp('BACKUP - LQR WITHOUT CONSTRAINTS');
    u = -param.K_LQR * (currentX - param.xTarget); 
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




formatSpec = 'MPC Controller ux=%f | uy=%f\n';
fprintf(formatSpec,u(1),u(2));
fprintf('------------------------------\n');
end % End of myMPController
