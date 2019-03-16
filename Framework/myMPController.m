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
if(param.selectController ~= 6)
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
% if(dist_final < param.closeToTarget || checkCloseTar == 1)
%    u = zeros(2,1);%stop do everything 
%    checkCloseTar = 1;
% end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Switch to LQR if close to the solution
             
% if(dist_final < param.closeToTarget || checkCloseTar == 1 && checkChangeRect == 1)
%     disp('USING LQR CLOSE TO THE FINAL SOLUTIO');
%     u = -param.K_LQR * (currentX - param.xTarget); 
%     checkCloseTar = 1;
% end
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
