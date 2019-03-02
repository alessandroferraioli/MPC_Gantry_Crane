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
