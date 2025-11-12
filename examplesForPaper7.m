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
runExample2_newtonIteration(2)
runExample2_newtonIteration(4)

% Fujimoto/Scherpen academic model (this one is nice)
runExample12_newtonIteration(2)
tic
runExample12_newtonIteration(4)
toc
tic
runExample12_newtonIteration(4)
toc

% Gray/Scherpen model
runExample13_newtonIteration(2)
runExample13_newtonIteration(4)

% Double pendulum gradient system
runExample14_newtonIteration(2)
runExample14_newtonIteration(4)

% Stable damped pendulum (this one is nice)
runExample31_newtonIteration(2)
runExample31_newtonIteration(4)

