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

persistent check
if(isempty(check))
   check = 0; 
end

dist = sqrt((x_hat(1)-param.xTarSigned)^2 + (x_hat(3)-param.yTarSigned)^2);
xStart = zeros(8,1);
if(check == 0)
    if(dist < param.epsilonTarget)
        %Second rect constr
        F = param.F2;
        J = param.J2;
        L = param.L2;
        EE = param.EE2;
        bb = param.bb2;
        newbb = bb - EE*kron(ones(param.N,1),uCurrentTarget);
        xStart(1) = param.xTarSigned;
        xStart(3) = param.yTarSigned;
        check = 1;
    else
        %First rect constr
        F = param.F;
        J = param.J;
        L = param.L;
        bb = param.bb;
        EE = param.EE;
        newbb = bb - EE*kron(ones(param.N,1),uCurrentTarget);
        xStart = param.xStart;
    end
else
    %Second rect constr
    F = param.F2;
    J = param.J2;
    L = param.L2;
    bb = param.bb2;
    EE = param.EE2;
    newbb = bb - EE*kron(ones(param.N,1),uCurrentTarget);
    xStart(1) = param.xTarSigned;
    xStart(3) = param.yTarSigned;
end


formatSpec = 'xHat =%f | yHat =%f | xT =%f | yT =%f | dtx =%f | dty =%f | utx =%f | uty =%f \n';
fprintf(formatSpec,currentX(1),currentX(3),r(1),r(3) ,dHat(1) ,dHat(2) , uCurrentTarget(1),  uCurrentTarget(2))

f =  (param.G * (currentX-xCurrentTarget)); %linear term must be a column vector
RHS = newbb+ L*xCurrentTarget+  J*xStart; %RHS of inequality
iA = false(size(bb));
[U,~,~]=mpcqpsolver(param.H,f,-F,-RHS,[],zeros(0,1),iA,opt);

u =U(1:param.m,:);

end % End of myMPController
