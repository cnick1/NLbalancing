% Examples for Paper 6 (by N. Corbin and A. Sarkar, in progress [1])
%
%  Description: This script runs the example and produces the data for 
%  reference [1]. The examples are:
%
%    - Example 3: Jeff's finite element Burgers' equation
%    - Example 6: a nonlinear (von Karman) Euler-Bernoulli beam with cable
%                 actuation that provides state dependent inputs. The drift 
%                 and input dynamics are up-to cubic. The model can be made
%                 arbitrarily large.
%    - Example 7: a popular 3D aircraft stall model from [2]. 
%    - Example 12: a cooked up model to demonstrate nonlinear balancing 
%                  transformations, originally from Fujimoto 2001 [3] but
%                  used several times in the literature.
%    - Example 13: a cooked up model to demonstrate nonlinear balancing 
%                  transformations, originally from Gray 2001 [4].
%    - Example 14: double pendulum model from Fujimoto 2008 [5].
%    - Example 15: a simple mass-spring-damper system with cubic stiffness;
%                  can be made arbitrarily large.
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
%             [5] K. Fujimoto and D. Tsubakino, “Computation of nonlinear 
%             balanced realization and model reduction based on Taylor 
%             series expansion,” Systems & Control Letters, vol. 57, no. 4,
%             pp. 283–289, Apr. 2008, doi: 10.1016/j.sysconle.2007.08.015
% 

close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

% exportData = true;

%% Example 3: Burgers Equation
% runExample3_outputDiagonalization();


%% Example 6: FEM Nonlinear Beam
% runExample6_outputDiagonalization();


%% Example 7: Garrard 3D Airplane Stall, check diagonalization
% [sigmaSquared1] = runExample7_outputDiagonalization();
% [sigmaSquared2] = runExample7_outputDiagonalization();
% If you modify outputDiag, you can check the effect of the different
% solutions to the transformations
% norm(sigmaSquared1 - sigmaSquared2) 

%% Example 12: Fujimoto 2D cooked up model for NL Balancing
% runExample12();

%% Example 13: Gray 2D cooked up model for NL Balancing
% runExample13();

%% Example 14b: 2D double pendulum model gradient system
% runExample14_gradient();

%% Example 14c: Krener's 6D triple pendulum model
% runExample14_krener()

%% Example 14: 4D double pendulum model
clc
runExample14();

%% Example 15: Coupled nonlinear mass-spring-damper model
% runExample15(32); % change model size n here



