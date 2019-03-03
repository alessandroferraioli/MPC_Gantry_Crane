function r = myTargetGenerator(x_hat, param)
%% Modify the following function for your target generation
dHat = x_hat(9:10);
xHat = x_hat(1:8);

H = zeros(10,10);
H(9,9) = 1; %min on ux
H(10,10) = 1; %min on uy
f = zeros(10,1);
r = zeros(10,1);
persistent changed

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    
    
    
    Aeq = [[eye(8) - param.A , -param.B];[param.Mx , zeros(2,2)]];
    beq = [param.Bd * dHat ; [xTar ; yTar] - param.Md*dHat ];
    
    Aineq = [-eye(10) ; eye(10)];
    low = [-xTar+param.eps_t 0 -yTar+param.eps_t 0 param.angCond 0 param.angCond 0 -param.ul']';
    high = [xTar+param.eps_t param.eps_r yTar+param.eps_t param.eps_r param.angCond param.eps_r param.angCond param.eps_r param.uh']';
    bineq = [low+[param.Cd*dHat ; 0 ; 0] ; high-[param.Cd*dHat ; 0 ; 0] ];
    
    [r,~,flag] = quadprog(H,f,Aineq,bineq,Aeq,beq,[],[],[]);
    
    
    if(flag == -2)
        [r,~,flag] = quadprog(H,f,Aineq,bineq,[],[],[],[],[]);
        if(flag == -2)
            Aineq = [zeros(10,8) , [zeros(8,2) ; -eye(2)] ; [zeros(10,8) , [zeros(8,2) ; eye(2)]]];
            bineq = [zeros(8,1) ; -param.ul ; zeros(8,1) ; param.uh];
            [r,~,~] = quadprog(H,f,Aineq,bineq,[],[],[],[],[],[]);
        end
        
    end
    
end
end % End of myTargetGenerator



