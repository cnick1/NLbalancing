function runExample32_balancingTransformation(degree,lim)
%runExample32_balancingTransformation Runs the 2D pendulum example to visualize the nonlinear balancing transformations.
%
%   Usage:  runExample32_balancingTransformation(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description: This simple 3D examples captures the key idea in the
%   model reduction problem: the presence of a subsystem that in some sense
%   contributes little (perhaps is decoupled) to the overall dynamics, yet
%   drives interactions that cannot directly be eliminated.
%       Consider the linear model from [1,2]:
%           ẋ₁ = −x₁ + 100 x₃ + u,
%           ẋ₂ = −2 x₂ + 100 x₃ + u,
%           ẋ₃ = −5 x₃ + u,
%            y = x₁ + x₂ + x₃,
%   The third state component is decoupled and decays quickly, so we
%   intuitively expect that we should be able to approximate this model
%   with a 2D model. However, x₃ strongly drives the states x₁ and x₂. This
%   in some sense directly demonstrates the need for balancing: the state
%   contributes little to observability (since it decays quickly and
%   contributes little to the output) but contributes significantly to
%   controllability (since it drives the other states).
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). The mapping from the z̄ coordinates
%   to the x coordinates is visualized by forming a grid in the z̄ coordinates
%   and mapping that grid to the x coordinates.
%
%   References: [3] P. Holmes, J. L. Lumley, G. Berkooz, and C. W. Rowley,
%                   Turbulence, coherent structures, dynamical systems and
%                   symmetry. Cambridge University Press, 2012. doi:
%                   10.1017/cbo9780511919701.
%
%   Part of the NLbalancing repository.
%%
% close all;
set(groot,'defaultLineLineWidth',1,'defaultTextInterpreter','TeX')

fprintf('Running Example 32\n')

if nargin < 2
    lim = 3.25;
    if nargin < 1
        degree = 4;
    end
end

%% Get system dynamics
[f, g, h] = getSystem32(transform=true);

%% Compute balanced realization
[fbal,gbal,hbal,Tbal] = getBalancedRealization(f,g,h,eta=0,degree=degree-1);
TbalInv = transformationInverse(Tbal);

%% Compute input-normal/output-diagonal realization
[v] = approxPastEnergy(f, g, h, 0, degree);
[w] = approxFutureEnergy(f, g, h, 0, degree);
[~, TinOd] = inputNormalOutputDiagonalTransformation(v, w, degree-1);
[finOd,ginOd,hinOd] = transformDynamics(f,g,h,TinOd);
[vbal, wbal] = transformEnergyFunctions(v,w,Tbal);
[vinOd, winOd] = transformEnergyFunctions(v,w,TinOd);

fprintf("\n  - FOM dynamics:\n\n")
dispKronPoly(f)

fprintf("\n  - Balanced dynamics:\n\n")
dispKronPoly(fbal,degree=degree)

fprintf("\n  - Energy Functions:\n\n")
dispKronPoly(v,n=3),dispKronPoly(w,n=3)

fprintf("\n  - Input-normal/output-diagonal energy Functions:\n\n")
dispKronPoly(vinOd,n=3),dispKronPoly(winOd,n=3)

fprintf("\n  - Balanced energy Functions:\n\n")
dispKronPoly(vbal,n=3),dispKronPoly(wbal,n=3)


end

