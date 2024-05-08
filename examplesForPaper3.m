% Examples for Paper 3 (by N. Corbin, submitted to IEEE TAC)
%
%  Description: This script runs all of the examples and produces all of
%  the data for reference [1], submitted to IEEE TAC. The examples are:
%
%    - Example 1: a 1-DOF quadratic-polynomial model for which we can
%      analytically plot the energy functions
%    - Example 2: a 2-DOF quadratic-bilinear model based on
%      Scherpen/Kawano [2]
%    - Example 6: a nonlinear (von Karman) Euler-Bernoulli beam with cable
%      actuation that provides state dependent inputs. The drift and input
%      dynamics are up-to cubic. The model can be made arbitrarily large.
%
%  Authors:   Nick Corbin, UCSD
% 
%  Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%             energy functions for polynomial control-affine systems,” 2023.
%             [2] Y. Kawano and J. M. A. Scherpen, “Model reduction by
%             differential balancing based on nonlinear hankel operators,”
%             IEEE Transactions on Automatic Control, vol. 62, no. 7,
%             pp. 3293–3308, Jul. 2017, doi: 10.1109/tac.2016.2628201.
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

exportData = true;
%%  Example 1: 1D ODE with analytical energy functions
% This generates Figure 1 (the energy functions)
runExample1(3, 3, exportData);

% This generates Figures 2 and 3 (the error plots)
runExample1_regionOfAccuracy(exportData);

%%  Example 2: 2D Energy Functions
% This generates the data for Figure 4 (the energy functions)
runExample2_energyFunctionPlots(exportData);

% This generates Figure 5 (the HJB residuals)
runExample2_regionOfAccuracy_res(exportData);

%% Example 6: Finite Element Beam Convergence
% This generates the data for Tables II and III, which are plotted in
% Figures 7, 8, and 9. (Note: this takes a few hours on a good modern
% laptop with 16GB RAM)
runExample6(4, 4, exportData, 0.1); % run this first to get the larger IC xb
runExample6(4, 4, exportData, 0.01); % run this first to get the smaller IC xa
