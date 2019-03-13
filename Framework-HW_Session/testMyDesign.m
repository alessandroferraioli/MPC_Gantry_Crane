clear variables
clc
close all


%% Load the parameters for the Simscape
load('Params_Simscape.mat');
load('SSmodelParams.mat');


%% Create the shape to test on
testShape = generateShape();


%% Extract the student functions
 extractFunctions(['FunctionTemplate.m'], 1);
 extractFunctions(['extrinsicFunctionCalls.m'], 1);
% 

%% Declare other simulation parameters
f = 20;
Ts = 1/f;

xZero = testShape.start(1,1);
yZero = testShape.start(1,2);



%% Call the setup function for the student
param = mySetup(testShape.c,...
                testShape.start,...
                testShape.target,...
                testShape.eps_r,...
                testShape.eps_t);


%% Save the data for the simulation 
save('workspace.mat');


%% Open the model
simModel = 'SimscapeCrane_ClosedLoop';
open(simModel); 


%% Import the data into Simulink
mws = get_param(simModel, 'modelworkspace');
mws.DataSource = 'MAT-File';
mws.FileName = 'workspace';
mws.reload();


%% Setup the simulation time
set_param(bdroot, 'StopTime', num2str(25) );


%% Update the controller blocks
updateScriptBlockContents( slroot, [simModel, '/MPController'], fileread('ext_MPC.m') );
updateScriptBlockContents( slroot, [simModel, '/State_Estimator'], fileread('ext_MSE.m') );
updateScriptBlockContents( slroot, [simModel, '/Target_Generator'], fileread('ext_MTG.m') );


%% Run the actual simulation
sim(simModel);


%% Visualize the Course
[~, ~] = analyzeCourse( [], GantryCraneOutput, GantryCraneInput, testShape, r, 0, 1 );


%% Visualize the performance
ul=[-1; -1];
uh=[1; 1];
cl=[0; 0];
ch=[xRange(2); yRange(2)];
xTarget = [testShape.target(1); 0; testShape.target(2); 0; 0; 0; 0; 0];
GantryResponsePlot(GantryCraneOutput.time,GantryCraneInput.signals.values,...
    GantryCraneOutput.signals.values,ul,uh,cl,ch,[1 3],xTarget,'Nonlinear simulation');