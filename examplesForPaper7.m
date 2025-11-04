% Examples for Paper 7 (balance-and-reduce)
%
%  Description: This script runs the example and produces the data for
%  reference [1]. The examples are:
%
%    - Example 23: 3D model from Sam Otto's papers with strong nonlinear
%                  coupling with a state that otherwise decays quickly
%    - Example 24: 2D model inspired by Example 23 to more easily visualize
%                  and understand the linear vs nonlinear transformations.
%
%  Authors:   Nick Corbin, UCSD
%
%  Reference: [1]
%

close all; clear; clc;
setKroneckerToolsPath
addpath('utils')
addpath('examples')

%% All balancing examples
runExample2_balancingTransformation(2)  % 2D Kawano QB model
runExample2_balancingTransformation(4)  % 2D Kawano QB model

runExample12_balancingTransformation(2) % 2D Fujimoto/Scherpen model
runExample12_balancingTransformation(4) % 2D Fujimoto/Scherpen model

runExample13_balancingTransformation(2) % 2D Gray/Scherpen model 
runExample13_balancingTransformation(4) % 2D Gray/Scherpen model 

runExample14_balancingTransformation(2) % 2D double pendulum gradient system 
runExample14_balancingTransformation(4) % 2D double pendulum gradient system 

runExample31_balancingTransformation(2) % 2D stable pendulum system 
runExample31_balancingTransformation(4) % 2D stable pendulum system 

%% OLD 
runExample11_balancingTransformation() % 2D inverted pendulum (REVISIT THIS WITH PPR EX 26!) ALSO STABLE PENDULUM! ALSO CHANGE THE ETA TO STABILIZE AND PLOT THE CLOSED-LOOP TRAJECTORIES

runExample15_balancingTransformation() % 4D double pendulum *** THIS IS THE FUJIMOTO/TSUBAKINO COMPARISON TO LOOK AT 
runExample23_balancingTransformation() % 3D Otto model
runExample24_balancingTransformation() % 2D similar to Otto model

runExample26_balancingTransformation(2)
runExample26_balancingTransformation(4)
runExample26_balancingTransformation(6)
runExample26_balancingTransformation(8)
runExample26_balancingTransformation(10)


%% Example 24: 2D models
% Before we look at a 3D model, let's look at a simple 2D model where we
% can see what the polynomial transformations are doing to our coordinate
% system. It seems to me that we really want some sort of bijective mapping
% between our original coordinates and our transformed coordinates. Our
% polynomial transformations are only bijective locally. What happens when
% they fail to be bijective?
runExample24(2)
runExample24(4)

% What does the whole space transformation look like?
runExample24_largerGrid(2); xlim([-.5 .5]);ylim([-.5 .5])
runExample24_largerGrid(4)

% What if we "zoom in" and look at a smaller grid
runExample24_smallerGrid(4)

%% Example 23: Sam Otto's 3D models
% There are a few things to compare here:
%   1) Even before we talk about truncation (reduction), what do the
%   polynomial transformations do regarding the solutions of the
%   differential equations governing the dynamics?

%   If we apply the transformation and invert the Jacobian at each point,
%   things seem to work "ok", although I have concerns about the
%   invertibility of the Jacobian:
runExample23(2,false,false,1) % Linear transformation

runExample23_jacobianInversion2(4,false,false,1) % applying the transformation at each step and inverting Jacobian at each step
%   The polynomial transformation *should* be able to be propagated through
%   the dynamics to get the explicit rhs before inverting the Jacobian.
%   However, it is quite expensive because then for example cubic dynamics
%   with a cubic transformation give degree 9 dynamics. I think it is
%   actually better to apply the transformation first then evaluate the
%   original dynamics. Currently it seems my f̃(z) = f(phi(z)) is not
%   evaluating properly, not sure if I did something wrong or if the
%   operation is just very poorly conditioned since I am blowing up to
%   large polynomials.
% runExample23_jacobianInversion(4,false,false,1) % Expansion approximation forward (should be exact), pointwise inversion backwards
%   The Jacobian inversion can also be approximated with a polynomial; this
%   would give explicit control-affine polynomial approximation to the
%   dynamics, but introduces some loss of accuracy (up to here in theory there was no truncation of series expansions, just roundoff error).
% runExample23(4,false,false,1) % Expansion approximation forward (should be exact), pointwise inversion backwards

%% Assuming we get the full transformation working properly, then
%% The question is:
%   2) When we introduce the truncation/reduction step, is the truncated
%   nonlinear transformation better than the truncated linear
%   transformation? There are several ways to compute the truncated reduced
%   order model (balance and reduce, etc.). The first way (true way) is to
%   compute the full transformation and then eliminate the last state in
%   the transformed model.

% These are the results from Sam's papers for LINEAR BT:
runExample23_balancingTransformation(2,false,true,1); ylim([0 50]) % Linear reduced transformation
runExample23_balancingTransformation(2,false,true,1/2); ylim([0 6]) % Linear reduced transformation


% What about nonlinear BT?
runExample23_balancingTransformation(4,false,true,1); ylim([0 50]) % Nonlinear reduced transformation
runExample23_balancingTransformation(4,false,true,1/2); ylim([0 6]) % Nonlinear reduced transformation

% runExample23(2,false,false,1) % Linear transformation
% runExample23_jacobianInversion2(2,false,false,1) % applying the transformation at each step and inverting Jacobian at each step

runExample23_balancingTransformation(2,false,true,1/50) % Linear reduced transformation
runExample23_balancingTransformation(4,false,true,1/50) % Nonlinear reduced transformation

runExample23_balancingTransformation(2,false,true,1/10) % Linear reduced transformation
runExample23_balancingTransformation(4,false,true,1/10) % Nonlinear reduced transformation

runExample23_balancingTransformation(2,false,true,1/2) % Linear reduced transformation
runExample23_balancingTransformation(4,false,true,1/2) % Nonlinear reduced transformation

runExample23_balancingTransformation(2,false,true,1) % Linear reduced transformation
runExample23_balancingTransformation(4,false,true,1) % Nonlinear reduced transformation






runExample23(2,true,true,1)

% Testing Jacobian expansion vs true pointwise Jacobian inversion
runExample23(4,false,false,1/10) % Expansion approximation
runExample23_jacobianInversion2(4,false,false,1/10) % "No approximations"
runExample23(2,false,false,1/10) % Linear transformation

% Still seems to work ok for larger initial condition


%% All balancing examples
runExample2_balancingTransformation()
runExample11_balancingTransformation()
runExample12_balancingTransformation()
runExample13_balancingTransformation()
runExample15_balancingTransformation()
runExample23_balancingTransformation()
runExample24_balancingTransformation()

%% On transient growth 
% It looks to me like our NLBT may not have the ability to do the thing
% with capturing transient growth? * Talking about the reduced model
runExample23_balancingTransformation(2,false,true,1/15) % Linear reduced transformation
runExample23_balancingTransformation(4,false,true,1/15) % Nonlinear reduced transformation

% Notice how the nonlinear ROM is not catching the transient growth, even
% though the linear projection does...
runExample23_balancingTransformation(2,false,true,1/10) % Linear reduced transformation
runExample23_balancingTransformation(4,false,true,1/10) % Nonlinear reduced transformation
