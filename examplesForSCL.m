%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Main script to run the examples for SCL %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Author: Nicholas Corbin, UCSD           %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Usage: Run the entire script to generate the examples for the SCL paper 
%         [1]. 
% 
%  Description: The examples used in reference [1] are:
%
%      Example 12: a cooked up model to demonstrate nonlinear balancing 
%                  transformations, originally from Fujimoto 2001 [2] but
%                  used several times in the literature.
%      Example 17: a simple mass-spring-damper system with cubic stiffness;
%                  can be made arbitrarily large.
%
%  Authors:   Nick Corbin, UCSD
%
%  Reference: [1]N. A. Corbin, A. Sarkar, J. M. A. Scherpen, and B. Kramer,
%             “Scalable computation of input-normal/output-diagonal
%             balanced realization for control-affine polynomial systems,”
%             Systems & Control Letters, vol. 204, p. 106178, Oct. 2025,
%             doi: 10.1016/j.sysconle.2025.106178.
%             [2] K. Fujimoto and J. M. A. Scherpen, “Model reduction
%             for nonlinear systems based on the differential
%             eigenstructure of Hankel operators,” in Proceedings of
%             the 40th IEEE Conference on Decision and Control (Cat.
%             No.01CH37228), IEEE, 2001. doi: 10.1109/cdc.2001.980322
% 
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

%% Example 12: Fujimoto 2D cooked up model for NL Balancing
runExample12();

%% Example 17: Coupled nonlinear mass-spring-damper model
% runExample17(32); % change model size n here

[ns,energyTimings,transformationTimings] = runExample17_timings(3);
[ns,energyTimings,transformationTimings] = runExample17_timings(4);
[ns,energyTimings,transformationTimings] = runExample17_timings(5);
[ns,energyTimings,transformationTimings] = runExample17_timings(6);

%%
runExample17_checkDiag(6)
runExample17_checkDiag(8)
runExample17_checkDiag(16)
runExample17_checkDiag(32)
runExample17_checkDiag(50)
runExample17_checkDiag(64)






