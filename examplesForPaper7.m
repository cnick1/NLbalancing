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
runExample2_newtonIteration(4)

% Fujimoto/Scherpen academic model (this one is nice)
runExample12_balancingTransformation(2)
runExample12_balancingTransformation(4)
runExample12_newtonIteration(4)


% Gray/Scherpen model
% This example is actually nice because the balanced realization is linear;
% we can probably work out what the original model was that they
% transformed using a nonlinear transformation. In fact, the balancing
% transformation computes the precise transformation we need ☺
runExample13_balancingTransformation(2)
runExample13_balancingTransformation(4)
runExample13_newtonIteration(4)

% Double pendulum gradient system
runExample14_balancingTransformation(2)
runExample14_balancingTransformation(4)
runExample14_newtonIteration(4)

% Stable damped pendulum (this one is nice)
runExample31_balancingTransformation(2)
runExample31_balancingTransformation(4)
runExample31_newtonIteration(4)

