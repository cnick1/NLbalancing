%  A script to run the examples in QB NLBT
%  The examples are the following: 
%
%       Example 5: a 1-DOF quadratic-bilinear model for which we can
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

%%  runExample5
% TODO: Fix Example 5
runExample5

%%  runExample7 produces the plots for Fig. 2.  
% %  For Table 1, please use energyFunctionValidation in the tests directory.
% [v,w] = runExample7(6,true,false,6);
%  figure
%  copyobj(allchild(1),3); 
%  figure
%  copyobj(allchild(2),4);
% [v,w] = runExample2(6,true,false,6);

%%  runExample8 produces the plots for Fig. 2.  
% [v,w] = runExample8(4,true);

