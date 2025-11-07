% Examples for Paper 7 (balance-and-reduce)
%
%  Description: This script runs the example and produces the data for
%  reference [1]. 
%
%  Authors:   Nick Corbin, UCSD
%
%  Reference: [1]
%

close all; clear; clc;
setKroneckerToolsPath
addpath('utils')
addpath('examples')

%% 2D Grid Transformation Examples
% Kawano QB model
runExample2_balancingTransformation(2)
runExample2_balancingTransformation(4)

% Fujimoto/Scherpen academic model (this one is nice)
runExample12_balancingTransformation(2)
runExample12_balancingTransformation(4)

% Gray/Scherpen model 
runExample13_balancingTransformation(2) 
runExample13_balancingTransformation(4) 

% Double pendulum gradient system 
runExample14_balancingTransformation(2) 
runExample14_balancingTransformation(4) 

% Stable damped pendulum (this one is nice)
runExample31_balancingTransformation(2) 
runExample31_balancingTransformation(4) 

