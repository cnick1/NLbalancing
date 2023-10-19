%  A script to run the examples in Nonlinear Balanced Truncation:
%  Part 1-Computing energy functions, Kramer, Gugercin, and Borggaard,
%  arXiv:2209.07645

setKroneckerToolsPath
addpath('examples')
addpath('utils')

%%  runExample1 produces data for Fig. 1.
runExample1

%%  runExample2 produces the plots for Fig. 2.  
%  For Table 1, please use runExample2_energyFunctionValidation
[v,w] = runExample2(6,true,false,6);

%%  runExample3 produces Tables II, III, and IV.
runExample3

%%  runExample4 produces Tables V and VI.
runExample4
