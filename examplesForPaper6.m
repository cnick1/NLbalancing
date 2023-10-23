% Examples for Paper 6 (by N. Corbin and A. Sarkar, in progress)
%
%  Description: This script runs the example and produces the data for 
%  reference [1]. The examples are:
%
%    - Example 8: a nonlinear heat equation from [2-4]. The dynamics are 
%      cubic. The model can be made arbitrarily large with finite elements.
%
%  Authors:   Nick Corbin, UCSD
%
%  Reference: [1] N. A. Corbin, A. Sarkar, and B. Kramer, “Computation of 
%             Output Diagonal Transformations for Nonlinear Balancing,” 
%             in progress, 2023.
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
addpath('ArijitExamples')
addpath('utils')

% exportData = true;

%% Checking input-normal/output-diagonal transformation
[v, w] = runExample7_outputDiagonalization();


