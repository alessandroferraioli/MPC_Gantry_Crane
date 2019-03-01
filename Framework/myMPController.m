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
