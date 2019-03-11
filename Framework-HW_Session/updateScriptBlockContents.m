function [ ] = updateScriptBlockContents( sl, blockPath, cont )
%UPDATEBLOCKCONTENTS Update script block with new script contents
%
% This function will update the contents inside a MATLAB script block
% inside a Simulink model. It will entirely replace the contents of the
% block.
%
% Usage:
%   updateScriptBlockContents( sl, blockPath, cont );
%
% Inputs:
%   sl        - Handle to the root Simulink diagram (found using slroot usually)
%   blockPath - The path inside the Simulink diagram to the block
%   cont      - A character vector containing the new contents for the block
%
% Example:
%   To update the block "controller" inside the model "clmodel" with the
%   contents of the file "controller.m" run the code:
%       updateScriptBlockContents( slroot, 'clmodel/controller', fileread('controller.m') );
%   This assumes that controller.m is in the current directory, and that
%   the Simulink model is already open.
% 

% Find the handle to the block to update
hand = find(sl, '-isa', 'Stateflow.EMChart', 'Path', blockPath);

% Error test the block handles
if ( isempty(hand) )
    error( ['Unable to get block handle for ', blockPath] );
end

% Read the controller, state estimator, and target generator into Simulink
hand.Script = cont;

% Get rid of the handle
clear hand;

end

