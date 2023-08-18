%  A script to run the examples in ACC
%  The examples are the following:
%
%       Example 8: a nonlinear heat equation. The dynamics are cubic. The model can be made arbitrarily large.
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

% exportData = false;
exportData = true;

%% Example 8: Finite Element Heat Nonlinear Equation
[w] = runExample8(exportData,1e-5);
