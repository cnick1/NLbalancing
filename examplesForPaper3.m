%  A script to run the examples in QB NLBT
%  The examples are the following: 
%
%       Example 1: a 1-DOF quadratic-bilinear model for which we can
%       analytically plot the energy functions 
%       Example 6: a 3-DOF quadratic-bilinear unicycle model 
%       Example 7: a 2-DOF quadratic-bilinear model based on
%       Scherpen/Kawano
%       Example 8: a nonlinear (von Karman) Euler-Bernoulli beam with cable
%       actuation that provides state dependent inputs. The dynamics and
%       input are up-to cubic. The model can be made arbitrarily large.
%       
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

exportData = true;
%%  Example 1: 1D ODE with analytical energy functions
% runExample1(1);
runExample1(2,2,exportData);
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
% plotEnergyFunctionsTest

%%  Example 2: 2D Energy Functions 
% [v,w] = runExample2(8,true,false,6,1,1,false,false);
% 
% % %  For Table 1, please use energyFunctionValidation in the tests directory.
% [v,w] = runExample2(6,true,false,6,1,1,false,false);
% [v,w] = runExample2(6,true,false,6,2,2,false,false);
% [v,w] = runExample2(6,true,false,6,3,3,false,false);
% 
% 
% % Kawano model
% [v,w] = runExample2(6,true,false,6,1,1,false,true);
[v,w] = runExample2(6,true,false,6,2,2,exportData,true);
% [v,w] = runExample2(6,true,false,6,3,3,false,true);

%% Example 2: Singular Value Functions
runExample2_singularValueFunctions(10, true, true, 2, 1, 1, false, false);

%%  runExample8 produces the plots for Fig. 2.  
% [v,w] = runExample8(4,true);

