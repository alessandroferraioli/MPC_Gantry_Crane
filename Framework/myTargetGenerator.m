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



