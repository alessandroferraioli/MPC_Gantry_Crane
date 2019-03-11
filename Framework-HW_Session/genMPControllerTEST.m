function [u,status,iA1] = genMPControllerTEST(H,G,F,bb,J,L,x,xTarget,m,iA)
% H       - quadratic term in the cost function (Linv if using mpcqpsolver).
% G       - matrix in the linear term in the cost function.
% F       - LHS of the inequalities term.
% bb      - RHS of the inequalities term.
% J       - RHS of inequalities term.
% L       - RHS of inequalities term. Use this to shift constraints around target point
% x       - current state
% xTarget - target equilibrium point.
% m       - Number of inputs.
% iA      - active inequalities, see doc mpcqpsolver
%
% u is the first input vector u_0 in the sequence [u_0; u_1; ...; u_{N-1}]; 
% In other words, u is the value of the receding horizon control law
% evaluated with the current state x0 and target xTarget

% Please read the documentation on mpcqpsolver to understand how it is
% suppose to be used. Use iA and iA1 to pass the active inequality vector 

opt = mpcqpsolverOptions;
opt.IntegrityChecks = false;%% for code generation
opt.FeasibilityTol = 1e-3;
opt.DataType = 'double';
%% your code starts here
f =  (G * (x-xTarget)); %linear term must be a column vector
RHS = bb+ L*xTarget+  J*x; %RHS of inequality 

%U has size m*N = 2*20 = 40
[U,status,iA1]=mpcqpsolver(H,f,-F,-RHS,[],zeros(0,1),iA,opt);
%% your remaining code here
u =U(1:m,:) ;
end

% Remember to copy and paste your code above to the Simulink model.