%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Balance / Reduce Paper Examples         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Author: Nicholas Corbin, UCSD           %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Usage: Run the entire script to generate the examples
%
%  Description: The leading examples I plan on using are
%     Ex. 6  - Cubic FEM beam
%     Ex. 12 - 2D Scherpen/Fujimoto cooked up example
%     Ex. 15 - 4D double pendulum
%     Ex. 31 - 2D stable pendulum (Newman)
%     Ex. 32 - 3D toy example w/ nonlinear transformation (Holmes)
%
%  Authors:   Nick Corbin, UCSD
%
%  Reference: [1]
%

close all; clear; clc;
setKroneckerToolsPath
addpath('utils')
addpath('examples')

%% 2D Fujimoto/Scherpen OR 2D Gray/Scherpen academic model 
% The purpose of this example is to illustrate the individual steps in the
% balancing transformation process. Since the dynamics are simple enough, we can
% display the equations for each of the steps in the process: the original
% analytic dynamics, the polynomial approximation to the dynamics, the nonlinear
% energy functions in the original coordinates, the transformations and the
% squared singular value functions, and ultimately the balanced realization.
% This example is interesting because, while it has been used many times in the
% literature and perhaps one could have suspected (I did not), the balanced
% realization is linear! So this example was constructed by taking a linear
% example and artifically applying a nonlinear transformation to it. This is an
% interesting thing to share with the reader, and the example otherwise also
% serves as a nice vehicle to show the grid transformations to illustrate the
% rotation/stretching of the nonlinear transformation using the polynomial
% approach vs the Newton iteration approach.
runExample12_nonlinearBalancing()
runExample13_nonlinearBalancing()

%% 2D Newman stable pendulum
% Continuing with the story, this example becomes one that is slightly more
% interesting because it is actually physical; the purpose here is to address
% the bijectivity of the transformation BEFORE truncation, so we can discuss
% that topic. We can reproduce the results of Newman, i.e. plotting the
% controllability and observability energies and their residuals, and plotting
% the system's response in the original coordinates vs the transformed
% coordinates. One of the purposes of this example is to highlight the low bar
% that has existed for other nonlinear balancing works: our method handles this
% example with ease, but in Newman's work this example already exhibited
% noticeable errors before any truncation of the model; thus, the proposed
% approach at least meets the wishlist item of having a locally bijective
% transformation before any truncation takes place. However, we can also
% highlight that this result from Newman's paper is in the region where the
% transformation is basically linear, so then perhaps we can have a result
% further away showing the benefit of the linear transformation: if it is
% locally bijective, it is globally bijective, whereas the nonlinear
% approach-even before truncation-introduces a locality result!
runExample31_nonlinearBalancing()

%% Holmes 3D model w/ nonlinear transformation
% Until now, none of the examples have considered truncation; this example
% serves to illustrate the projection of the ROM onto a manifold. I can either
% state up-front or at the end that the FOM here is obtained by transforming
% Holmes' linear example with a nonlinear transformation; it will just depend on
% how the story flows. The main result will be the manifold figure with the
% trajectories, and a plot of the true output, the linear ROM output, and the
% nonlinear ROM output with corresponding L2 error metrics. 
runExample32_balancedReduction(2)
runExample32_balancedReduction(4)
runExample32_balancedReduction(4,'white-noise')

%% 4D Double Pendulum
runExample15_nonlinearBalancing()

%% Cubic FEM beam
runExample6_nonlinearBalancing()

