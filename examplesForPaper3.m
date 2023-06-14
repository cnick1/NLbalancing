%  A script to run the examples in QB NLBT
%  The examples are the following:
%
%       Example 1: a 1-DOF quadratic-bilinear model for which we can
%       analytically plot the energy functions
%       Example 2: a 2-DOF quadratic-bilinear model based on
%       Scherpen/Kawano
%       Example 6: a nonlinear (von Karman) Euler-Bernoulli beam with cable
%       actuation that provides state dependent inputs. The dynamics and
%       input are up-to cubic. The model can be made arbitrarily large.
%       Example 5: a 3-DOF quadratic-bilinear unicycle model
%
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

% exportData = false;
exportData = true;
kawanoModel = true;
%%  Example 1: 1D ODE with analytical energy functions
% runExample1(1);
runExample1(3,3,exportData);
% runExample1(3);

% Compare all 4 cases
% a = 4; b = 5;
% runExample1(1,1);
% close; xlim([-a, a]); ylim([0, b]); grid on
% title('quadratic, quadratic assumed')
% runExample1(2,1); close; xlim([-a, a]); ylim([0, b]); grid on
% title('QB system, quadratic assumed')
% runExample1(2,2); close; xlim([-a, a]); ylim([0, b]); grid on
% title('QB system, QB assumed')
% runExample1(3,3); close; xlim([-a, a]); ylim([0, b]); grid on
% title('polynomial system, polynomial assumed ***')
%

runExample1_regionOfAccuracy(exportData);
runExample1_regionOfAccuracy_res(exportData);

%%  Example 2: 2D Energy Functions
% For Table 1, please use energyFunctionValidation in the tests directory.
% [v, w] = runExample2(6, true, false, 6, 1, 1, false, kawanoModel);
[v, w] = runExample2(4, true, false, 4, 2, 2, exportData, kawanoModel);
% [v, w] = runExample2(6, true, false, 6, 3, 3, false, kawanoModel);


% Example 2: Singular Value Functions
runExample2_singularValueFunctions(4, true, true, 1, 2, 2, exportData, kawanoModel);
runExample2_singularValueFunctions(4, true, true, 1, 2, 2, exportData, kawanoModel);
runExample2_singularValueFunctions(4, true, true, 3, 2, 2, exportData, kawanoModel);
close all

% Example 2: residuals


%% Example 6: Finite Element Beam Convergence
% [w] = runExample6(2,2,false);
[w] = runExample6(2,2,exportData,0.1);
[w] = runExample6(3,3,false,0.01);
[w] = runExample6(2,2,exportData,0.01);
