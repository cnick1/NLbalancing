%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main script to run the examples for ACC %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Author: Nicholas Corbin, UCSD           %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Usage: Run the entire script to generate the examples for the ACC 2024
%         Toronto paper [1]. 
% 
%  Description: The examples used in reference [1] are:
%
%      Example 8: a nonlinear heat equation from [2-4]. The dynamics are 
%      cubic. The model can be made arbitrarily large with finite elements.
%
%  Authors:   Nick Corbin, UCSD
%
%  Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%             energy functions for polynomial drift nonlinear systems,” in
%             2024 American Control Conference (ACC), IEEE, Jul. 2024, pp.
%             2506–2511. doi: 10.23919/acc60939.2024.10644363.
%             [2] M. Embree, “Unstable modes in projection-based
%             reduced-order models: how many can there be, and what do they
%             tell you?,” Systems & Control Letters, vol. 124, pp. 49–59,
%             Feb. 2019, doi: 10.1016/j.sysconle.2018.11.010
%             [3] J. Galkowski, “Nonlinear instability in a semiclassical
%             problem,” Communications in Mathematical Physics, vol. 316,
%             no. 3, pp. 705–722, Oct. 2012, doi: 10.1007/s00220-012-1598-5
%             [4] B. Sandstede and A. Scheel, “Basin boundaries and
%             bifurcations near convective instabilities: a case study,”
%             Journal of Differential Equations, vol. 208, no. 1, pp.
%             176–193, Jan. 2005, doi: 10.1016/j.jde.2004.02.016
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

% If set to false, the code will only run the examples that are small
% enough to run on 16GB ram; set to true to run all of the examples in the
% paper.
exportData = false; 

%% Example 8: Finite Element Heat Nonlinear Equation
runExample8(exportData, 5e-5);
