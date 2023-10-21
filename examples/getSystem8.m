function [f, g, h] = getSystem8(numElements, elementOrder)
%getSystem8  Generates a cubic finite element model system for testing
%            energy functions. The system is a finite element model for a
%            nonlinear heat equation, i.e. a reaction-diffusion equation.
%            The function returns a finite element model with either linear
%            or quadratic elements (elementOrder 1 or 2). The state space
%            dimension is n=numElements-1 for linear elements or
%            n=2*numElements-1 for quadratic elements, due to the fixed end
%            boundary conditions.
%
%   Usage:   [f,g,h] = getSystem8()
%         or [f,g,h] = getSystem8(numElements,elementOrder)
%
%   Inputs:
%       numElements    - number of elements to discretize the domain with
%                        (default = 4)
%       elementOrder   - select either:
%                          • 1 = linear elements (default)
%                          • 2 = quadratic elements (not fully implemented)
%
%   Outputs: f,g,h     - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output (generalizations
%                        containing A,B,C, and N)
%
%   Background: after finite element discretization, the finite element
%   equations for the reaction-diffusion problem can be written as
%
%     M ẋ + K(x) x = B(x) u,
%     y = C x.
%
%   The terms K(x) and B(x) can be approximated with Taylor series
%   expansions, which leads to the Kronecker product representation
%
%     M ẋ + K₁ x + K₂ (x ⊗ x) + K₃ (x ⊗ x ⊗ x) + ...
%       = B₀ u + B₁ (x ⊗ u) + B₂ (x ⊗ x ⊗ u) + ...,
%     y = C x.
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial drift nonlinear systems,” 2023.
%              [2] M. Embree, “Unstable modes in projection-based
%               reduced-order models: how many can there be, and what do they
%               tell you?,” Systems & Control Letters, vol. 124, pp. 49–59,
%               Feb. 2019, doi: 10.1016/j.sysconle.2018.11.010
%              [3] J. Galkowski, “Nonlinear instability in a semiclassical
%               problem,” Communications in Mathematical Physics, vol. 316,
%               no. 3, pp. 705–722, Oct. 2012, doi: 10.1007/s00220-012-1598-5
%              [4] B. Sandstede and A. Scheel, “Basin boundaries and
%               bifurcations near convective instabilities: a case study,”
%               Journal of Differential Equations, vol. 208, no. 1, pp.
%               176–193, Jan. 2005, doi: 10.1016/j.jde.2004.02.016
%
%   Part of the NLbalancing repository.
%%

vec = @(X) X(:);

if nargin < 1
    numElements = 4;
end

if nargin < 2
    elementOrder = 1;
end

if elementOrder == 2
    error("Not fully implemented!")
end

rodLength = 30;
alpha = -1/8;

%% Define geometry and properties

% Define element properties
numNodes = elementOrder * numElements + 1; % number of nodes
x = linspace(0, rodLength, numNodes); % node locations
elementLength = x(elementOrder + 1) - x(1); % element length

% Define DOF counts
DOFsPerNode = 1;
DOFsPerElement = (elementOrder + 1) * DOFsPerNode;
TotalDOFs = numNodes * DOFsPerNode;

%% Assemble linear global matrices (mass and stiffness)
% Define mass matrix for one element
if elementOrder == 1
    M1E = elementLength / 6 * ...
        [2, 1;
     1, 2];

    % Define stiffness matrix for one element
    K1E = 1 / elementLength * [1, -1;
                               -1, 1] ...
        + alpha * elementLength / 6 * [2, 1;
                                   1, 2] ...
        +1/2 * [-1 1;
            -1 1];
elseif elementOrder == 2
    M1E = elementLength / 30 * ...
        [4, 2, -1;
     2, 16, 2;
     -1, 2, 4];

    % Define stiffness matrix for one element
    K1E = 1 / (3 * elementLength) * [7, -8, 1;
                                     -8, 16, -8;
                                     1, -8, 7] ...
        - elementLength / 240 * [4, 2, -1;
                             2, 16, 2;
                             -1, 2, 4];
end
% Initialize and stack/assemble global matrix
M1G = sparse(TotalDOFs, TotalDOFs);
K1G = sparse(TotalDOFs, TotalDOFs);

if elementOrder == 1
    for ii = 1:numElements
        M1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) = M1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) + M1E;
        K1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) = K1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) + K1E;
    end
elseif elementOrder == 2
    for i = 1:numElements
        ii = 2 * i - 1;
        M1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) = M1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) + M1E;
        K1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) = K1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) + K1E;
    end
end

%% Assemble quadratic global matrix
% Form all zero quadratic global matrix directly
K2G = sparse(TotalDOFs, TotalDOFs ^ 2);

%% Assemble cubic global matrix
% Define stiffness matrix for one element
if elementOrder == 1
    K3E = -elementLength / 20 * [4, 3, 0, 2, 0, 0, 0, 1;
                                 1, 2, 0, 3, 0, 0, 0, 4];
elseif elementOrder == 2
    K3E = -elementLength / 1260 * ...
        [92, 96, -21, 0, 96, -24, 0, 0, 6, 0, 0, 0, 0, 32, -48, 0, 0, -12, 0, 0, 0, 0, 0, 0, 0, 0, -7;
     32, 96, -12, 0, 96, -96, 0, 0, -12, 0, 0, 0, 0, 512, 96, 0, 0, 96, 0, 0, 0, 0, 0, 0, 0, 0, 32;
     -7, -12, 6, 0, -48, -24, 0, 0, -21, 0, 0, 0, 0, 32, 96, 0, 0, 96, 0, 0, 0, 0, 0, 0, 0, 0, 92];
end

% Initialize and stack/assemble global matrix
K3G = sparse(TotalDOFs, TotalDOFs ^ 3);

if elementOrder == 1
    for i = 0:numElements - 1 % start from zero so you don't have to subtract 1 every time
        ii = i * DOFsPerNode + 1;
        idxs = (TotalDOFs ^ 2 + TotalDOFs + 1) * i * DOFsPerNode ... % Starting index shift depending on element iteration
            + vec(( ...
            vec(([1:DOFsPerElement] + [0:TotalDOFs:TotalDOFs * (DOFsPerElement - 1)]')')' ... % (basically the quadratic indices)
            + [0:TotalDOFs ^ 2:TotalDOFs ^ 2 * (DOFsPerElement - 1)]' ... % Add secondary skips into sequence (add row to column and then vec)
        )')';

        % "stack" element matrices into global matrix
        K3G(ii:(ii + DOFsPerElement - 1), idxs) = K3G(ii:(ii + DOFsPerElement - 1), idxs) + K3E;
    end
elseif elementOrder == 2
    for i = 1:numElements
        ii = 2 * i - 1;
        idxs = (TotalDOFs ^ 2 + TotalDOFs + 1) * (i - 1) * DOFsPerNode * 2 ... % Starting index shift depending on element iteration
            + vec(( ...
            vec(([1:DOFsPerElement] + [0:TotalDOFs:TotalDOFs * (DOFsPerElement - 1)]')')' ... % (basically the quadratic indices)
            + [0:TotalDOFs ^ 2:TotalDOFs ^ 2 * (DOFsPerElement - 1)]' ... % Add secondary skips into sequence (add row to column and then vec)
        )')';

        % "stack" element matrices into global matrix
        K3G(ii:(ii + DOFsPerElement - 1), idxs) = K3G(ii:(ii + DOFsPerElement - 1), idxs) + K3E;
    end
end

%% RHS
% TODO
% RB0 = sparse(TotalDOFs, 2);
%         RB1 = sparse(TotalDOFs, 2 * TotalDOFs);
%         RB2 = sparse(TotalDOFs, 2 * TotalDOFs ^ 2);
%         RB3 = sparse(TotalDOFs, 2 * TotalDOFs ^ 3);

%         RB0(TotalDOFs - 2, 2) = 1; % Force in x direction
%         RB0(TotalDOFs - 1, 1) = 1; % Force in y direction
%         RB0(TotalDOFs, :) = 0; % Moment in z direction

RB0 = generate_B_matrix(4, (TotalDOFs - 1) / 4);

%% Impose boundary conditions
fixedDOFs = [1, TotalDOFs]; % First and last nodes fixed
freeDOFs = setdiff(1:TotalDOFs, fixedDOFs);

% Reduced system method
K1G = K1G(freeDOFs, freeDOFs);
M1G = M1G(freeDOFs, freeDOFs);
RB0 = RB0(freeDOFs, :);

% K2G could clean up
fixedDOFsSquared = vec(fixedDOFs.' + (0:TotalDOFs:TotalDOFs ^ 2 - 1));
fixedDOFsSquared = unique([fixedDOFsSquared; vec((1:TotalDOFs).' + (fixedDOFs - 1) * TotalDOFs)]);

freeDOFsSquared = setdiff(1:TotalDOFs ^ 2, fixedDOFsSquared);
K2G = K2G(freeDOFs, freeDOFsSquared);

% K3G could clean up
fixedDOFsCubed = vec(vec((fixedDOFs - 1) * TotalDOFs + [1:TotalDOFs].') + (0:TotalDOFs ^ 2:TotalDOFs ^ 3 - 1));
fixedDOFsCubed = [fixedDOFsCubed; vec((1:TotalDOFs ^ 2).' + (fixedDOFs - 1) * TotalDOFs ^ 2)]; % Top rows zero
fixedDOFsCubed = [fixedDOFsCubed; vec(vec(fixedDOFs.' + (0:TotalDOFs:TotalDOFs ^ 2 - 1)) + (0:TotalDOFs ^ 2:TotalDOFs ^ 3 - 1))];
fixedDOFsCubed = sort(unique(fixedDOFsCubed));

freeDOFsCubed = setdiff(1:TotalDOFs ^ 3, fixedDOFsCubed);
K3G = K3G(freeDOFs, freeDOFsCubed);

%% Convert to state-space representation
n = length(M1G);

McholL = chol(M1G).'; % Use Cholesky factor for inverting rather than inv()

F1 = -McholL.' \ (McholL \ K1G);

G0 = McholL.' \ (McholL \ RB0);

C = RB0.';

% Construct N₂
F2 = -McholL.' \ (McholL \ K2G);

% Construct N₃
F3 = -McholL.' \ (McholL \ K3G);

%% Format outputs
f = {full(F1), F2, F3};
g = {full(G0)};
h = {full(C)};

A = full(F1);
B = full(G0);
C = full(C);
N = full(F2);

end

function B = generate_B_matrix(m, i)
if rem(i, 1) ~= 0
    error("Incorrect number of elements/inputs")
end
n = i * m + 1;
B = sparse(n, m);
for j = 0:m - 1
    B(i * j + 1:i * j + 1 + i, j + 1) = 1;
end

B = B ./ (i + 1);
end
