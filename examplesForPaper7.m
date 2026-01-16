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

%% 2D Fujimoto/Scherpen academic model
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
% This example is actually very interesting and tells an important part of the
% story. PRIOR to truncation, the linear transformation is superior because it
% retains all of the information. However, when it comes to truncation, even for
% trajectories close to the equilibrium point, there is a qualitative difference
% between the linear and nonlinear transformations that results in the nonlinear
% model under the linear truncated transformation behaving more like the
% linearized dynamics than the nonlinear dynamics, even though the dynamics are
% still nonlinear! So I was wrong previously when I thought that the comparison
% shown was the linearized ROM vs the nonlinear ROM.
runExample15_balancedReduction(1,1,2,1)
runExample15_balancedReduction(1,3,2,1)
runExample15_balancedReduction(3,3,2,1)

% Illustration that error is lower locally but higher far away 
degrees = [1:2:5]; U0s = 1:.25:9;
errors = zeros(length(degrees),length(U0s));
for i = 1:length(degrees)
    for j = 1:length(U0s)
        errors(i,j) = runExample15_balancedReduction(degrees(i),5,2,U0s(j));
    end
end
figure; hold on;
for i = 1:length(degrees)
    plot(U0s,errors(i,:),'DisplayName',sprintf('Degree %i transformation',degrees(i)-1))
end
legend


%% Cubic FEM beam
% The story will wrap up centered around that plot of error vs distance to
% equilibrium: basically for all problems, that plot is qualitatively the
% behavior to expect. For problems like the Holmes example with a polynomial
% nonlinear transformation, we are able to move that intersection point far off
% to infinity, but in general there is no reason to expect that to be the case,
% and instead it will likely be near the equilibrium point of interest.
% *** I should make that plot for the Holmes example too! For "distance from
% equilibrium", what happens if I choose initial conditions all ON the balanced
% manifold, vs. not?
% ** The balanced manifold is not an invariant manifold! That's where the white
% noise excitation is a better indicator: certainly there are initial conditions
% for which the balanced ROM may perform poorly, and likewise where the linear
% ROM may do better; the issue is whether those initial conditions are
% reachable/observable! That's why, in the Holmes case, choosing an initial
% condition from the balanced manifold works well: it is
% controllable/observable. However, in general, the balanced manifold
% approximation we compute is not going to be exact, as it was for the Holmes
% example. Thus, for MSD, beam, or heat eq., can I analytically pick initial
% conditions from the balanced manifold? The trade-off is one of reducibility
% vs. numerical conditioning: the beam would be a nice example, but it yields
% poor numerical conditioning.
%
% In a twist, I've been able to get the beam example to work somewhat
% satisfactorily. The story is still the same: there is a trade-off between
% reducibility and numerical conditioning. Therefore, I was unable to use
% the model with just tip measurements and control. However, if I add
% control to each note and measurements to each node, the model does become
% better conditioned and is still moderately reducible. Curiously, I am
% seeing the behavior where the horizontal measurement is zero for the
% nonlinear model with linear transformation AFTER reduction... so the
% nonlinear output information is in certain states in the linear
% transformation and different states in the nonlinear transformation, such
% that that information is truncated in one case but not the other. THAT is
% VERY interesting. 
runExample6_nonlinearBalancing()

