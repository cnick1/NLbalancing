% Examples for Paper 6 (by N. Corbin and A. Sarkar, in progress [1])
%
%  Description: This script runs the example and produces the data for 
%  reference [1]. The examples are:
%
%    - Example 7: a popular 3D aircraft stall model from [2]. 
%    - Example 12: a cooked up model to demonstrate nonlinear balancing 
%      transformations, originally from Fujimoto 2001 [3] but used several 
%      times in the literature. 
%    - Example 13: a cooked up model to demonstrate nonlinear balancing 
%      transformations, originally from Gray 2001 [4].
%
%  Authors:   Nick Corbin, UCSD
%
%  Reference: [1] N. A. Corbin, A. Sarkar, and B. Kramer, “Computation of 
%             Output Diagonal Transformations for Nonlinear Balancing,” 
%             in progress, 2023.
%             [2] W. L. Garrard and J. M. Jordan, “Design of nonlinear
%             automatic flight control systems,” Automatica, vol. 13, no. 5,
%             pp. 497–505, Sep. 1977, doi: 10.1016/0005-1098(77)90070-x
%             [3] K. Fujimoto and J. M. A. Scherpen, “Model reduction
%             for nonlinear systems based on the differential
%             eigenstructure of Hankel operators,” in Proceedings of
%             the 40th IEEE Conference on Decision and Control (Cat.
%             No.01CH37228), IEEE, 2001. doi: 10.1109/cdc.2001.980322
%             [4] W. S. Gray and J. M. A. Scherpen, “On the nonuniqueness
%             of singular value functions and balanced nonlinear
%             realizations,” Systems & Control Letters, vol. 44, no. 3,
%             pp. 219–232, Oct. 2001, doi: 10.1016/s0167-6911(01)00144-x
% 

close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

% exportData = true;

%% Example 7: Garrard 3D Airplane Stall, check diagonalization
% [sigmaSquared] = runExample7_outputDiagonalization();


%% Example 12: Fujimoto 2D cooked up model for NL Balancing
runExample12()

%% Example 13: Gray 2D cooked up model for NL Balancing
runExample13()



