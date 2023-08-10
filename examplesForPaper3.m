%  A script to run the examples in Polynomial NLBT
%  The examples are the following:
%
%       Example 1: a 1-DOF quadratic-polynomial model for which we can
%       analytically plot the energy functions
%       Example 2: a 2-DOF quadratic-bilinear model based on
%       Scherpen/Kawano
%       Example 6: a nonlinear (von Karman) Euler-Bernoulli beam with cable
%       actuation that provides state dependent inputs. The drift and input
%       dynamics are up-to cubic. The model can be made arbitrarily large.
%
%   Reference: Submitted to IEEE TAC. 
%   "Scalable Computation of ℋ∞ Energy Functions for Polynomial 
%    Control-Affine Systems", N. Corbin and B. Kramer
%    arXiv:
% 

close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

exportData = true;
%%  Example 1: 1D ODE with analytical energy functions
% This generates Figure 1
runExample1(3,3,exportData); 

% This generates Figures 2 and 3
runExample1_regionOfAccuracy(exportData);

%%  Example 2: 2D Energy Functions
% This generates the data for Figure 4
runExample2(4, true, false, 4, 2, 2, exportData, true); 

% This generates Figure 5
runExample2_regionOfAccuracy_res(exportData);

%% Example 6: Finite Element Beam Convergence
% This generates the data for Tables II and III, which are plotted in
% Figures 7, 8, and 9. (Note: this takes a few hours on a good modern  
% laptop with 16GB RAM)
runExample6(4,4,exportData,0.1); % run this first to get the larger IC xb
runExample6(4,4,exportData,0.01); % run this first to get the smaller IC xa
