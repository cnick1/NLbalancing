%  A script to run the examples in QB NLBT
close all; clear; clc;
setKroneckerToolsPath
addpath('examples')
addpath('utils')

%%  runExample5
% runExample5


%%  runExample7 produces the plots for Fig. 2.  
% %  For Table 1, please use energyFunctionValidation in the tests directory.
[v,w] = runExample7(6,true,false,6);
%  figure
%  copyobj(allchild(1),3); 
%  figure
%  copyobj(allchild(2),4);
% [v,w] = runExample2(6,true,false,6);

vTest = load('testV.mat'); vTest = vTest.v;

%%  runExample8 produces the plots for Fig. 2.  
% [v,w] = runExample8(6,true,false,6);

