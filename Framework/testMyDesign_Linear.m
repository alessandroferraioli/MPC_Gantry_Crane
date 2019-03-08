clear variables;
close all;
clc;



%% Load the parameters for the Simscape
load('Params_Simscape.mat');
load('SSmodelParams.mat');


%% Create the shape to test on
testShape = generateShape();


%% Extract the student functions
extractFunctions(['FunctionTemplate.m'], 1);


%% Declare other simulation parameters
f = 20;
Ts = 1/f;



xZero = testShape.start(1,1);
yZero = testShape.start(1,2);

%% Set the simulation time
T = 20;
N=ceil(T/Ts);

%% Call the setup function for the student
param = mySetup(testShape.c,...
                testShape.start,...
                testShape.target,...
                testShape.eps_r,...
                testShape.eps_t);


%% Create the dynamics matrices
[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);
sysd=ss(A,B,C,0,Ts);




%% Run the actual linear simulation
% create waiting bar
hw=waitbar(0,'Please wait...');
warning('on');

% Initial conditions
x=[xZero; 0; yZero; 0; 0; 0; 0; 0];
y = x;
u = [0; 0];
u_all = [];


% Variable to hold optimization time
numSamples = T/Ts;
allContTime = [];

t=0:Ts:T;

%%
% Iterate for the simulation time
for t=0:Ts:T
    waitbar(t/T,hw,'Please wait...');
    tic;

    % Call the state estimator
    x_hat = myStateEstimator(u, y, param);

    % Call the target generator
    ref = myTargetGenerator(x_hat, param);
    
    % Call the controller function
    u = myMPController(ref, x_hat, param);
    
    contTime=toc;

    % Simulate
    [y, tt, xx] = lsim(sysd, [u';0 0], [0 Ts], x(:,end));

    % Keep the state variables around for analysis
    x=[x xx(end,:)'];
    u_all = [u_all u];
    
    % Save only the most recent output
    y = y(end,:)';
    
    % Save the computation time for analysis
    allContTime=[allContTime; contTime];
end

close(hw);
t=0:Ts:t;
x=x(:,1:length(t))';


%% Analyze the controller runtime
figure('Name','Controller Runtime');
plot(t, allContTime);
xlabel('Simulation time [s]')
ylabel('Controller runtime [s]')


%% Visualize the Course
GantryCraneOutput.time = t;
GantryCraneOutput.signals.values = x;
GantryCraneInput.signals.values = u_all';

[~, ~] = analyzeCourse( [], GantryCraneOutput, GantryCraneInput, testShape, r, 0, 1 );

%% Visualize the performance
ul=[-1; -1];
uh=[1; 1];
cl=[0; 0];
ch=[xRange(2); yRange(2)];
xTarget = [testShape.target(1); 0; testShape.target(2); 0; 0; 0; 0; 0];
GantryResponsePlot(GantryCraneOutput.time,GantryCraneInput.signals.values,...
    GantryCraneOutput.signals.values,ul,uh,cl,ch,[1 3],xTarget,'Linear simulation');
