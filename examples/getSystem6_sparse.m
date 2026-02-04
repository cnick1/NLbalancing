function [f, g, h] = getSystem6_sparse(numElements)
%getSystem6_sparse  Generates a cubic finite element beam model. The system models a
% nonlinear (due to von Karman strains) Euler-Bernoulli Beam with numElements
% elements, returning a state-space system with 6*numElements degrees of freedom
% (because each element has two nodes, each with 6 degrees of freedom (3
% position, 3 velocity), and the first node is fixed).
%
%   Usage:   [f,g,h] = getSystem6_sparse(numElements,actuatorConfig,rotaryInertia)
%
%   Inputs:
%       numElements    - number of elements to discretize the beam with
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients for
%                        the drift, input, and output in generalized form
%
%   Description: After discretization, the finite element equations for the beam
%   can be written as
%
%           M q̈ + D q̇ + K(q) q = B(q) u,
%           y = C₁ q̇ + C₂ q.
%
%   The terms K(q) and B(q) can be approximated with Taylor series
%   expansions, which leads to the Kronecker product representation
%
%     M q̈ + D q̇ + K₁ q + K₂ (q ⊗ q) + K₃ (q ⊗ q ⊗ q) + ...
%       = B₀ u + B₁ (q ⊗ u) + B₂ (q ⊗ q ⊗ u) + ...,
%     y = C₁ q̇ + C₂ q.
%
%   Defining the state vector x = [q q̇]^T, we can convert to a first-order
%   nonlinear state-space
%
%     ẋ = A x + N₂ (x ⊗ x) + N₃ (x ⊗ x ⊗ x) ...
%          + G₀ u + G₁ (x ⊗ u) + G₂ (x ⊗ x ⊗ u) + G₃ (x ⊗ x ⊗ x ⊗ u)
%     y  = C₁ q̇ + C₂ q.
%
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗∞
%               energy functions for polynomial control-affine systems,"
%               IEEE Transactions on Automatic Control, pp. 1–13, 2024,
%               doi: 10.1109/tac.2024.3494472
%
%   Part of the NLbalancing repository.
%%
arguments
numElements = 3
end
vec = @(X) X(:);
%% Define beam geometry and properties
BeamLength = 1; % length of beam
ElasticModulus = 210e9; % Young's modulus
CrossSecArea = 5e0; % cross-sectional area
% MomOfInertia = pi^2/16/1.8751040^2*CrossSecArea;
MomOfInertia = 1e-2; % moment of inertia .5e-2
density = 8000; % density
delta = 0.1; % Cable attachment distance from centerline of beam if actuatorConfiguration = 2
d1 = 0.0001; d2 = 0.0001;

% Define element properties
numNodes = numElements + 1; % number of nodes
nel = numElements;

% Define DOF counts
DOFsPerNode = 3;
DOFsPerElement = 2 * DOFsPerNode;
TotalDOFs = numNodes * DOFsPerNode;


%% Generate a mesh on a linear domain
% Get connectivity matrix and node locations
[mconn,x] = generateMesh(BeamLength,numElements);
a = x(2) - x(1); % element length

%% Assemble linear global matrices (mass and stiffness)
% Define mass matrix for one element
M1E = density * CrossSecArea * a / 420 * ...
      [140, 0, 0, 70, 0, 0;
      0, 156, 22 * a, 0, 54, -13 * a;
      0, 22 * a, 4 * a ^ 2, 0, 13 * a, -3 * a ^ 2;
      70, 0, 0, 140, 0, 0;
      0, 54, 13 * a, 0, 156, -22 * a;
      0, -13 * a, -3 * a ^ 2, 0, -22 * a, 4 * a ^ 2];

% Define stiffness matrix for one element
K1E = [ElasticModulus * CrossSecArea / a, 0, 0, -ElasticModulus * CrossSecArea / a, 0, 0;
      0, 12 * ElasticModulus * MomOfInertia / a ^ 3, 6 * ElasticModulus * MomOfInertia / a ^ 2, 0, -12 * ElasticModulus * MomOfInertia / a ^ 3, 6 * ElasticModulus * MomOfInertia / a ^ 2;
      0, 6 * ElasticModulus * MomOfInertia / a ^ 2, 4 * ElasticModulus * MomOfInertia / a, 0, -6 * ElasticModulus * MomOfInertia / a ^ 2, 2 * ElasticModulus * MomOfInertia / a;
      -ElasticModulus * CrossSecArea / a, 0, 0, ElasticModulus * CrossSecArea / a, 0, 0;
      0, -12 * ElasticModulus * MomOfInertia / a ^ 3, -6 * ElasticModulus * MomOfInertia / a ^ 2, 0, 12 * ElasticModulus * MomOfInertia / a ^ 3, -6 * ElasticModulus * MomOfInertia / a ^ 2;
      0, 6 * ElasticModulus * MomOfInertia / a ^ 2, 2 * ElasticModulus * MomOfInertia / a, 0, -6 * ElasticModulus * MomOfInertia / a ^ 2, 4 * ElasticModulus * MomOfInertia / a];

% Define stiffness matrix for one element
K2E = sparse([[0, 0, 0, 0, 0, 0;
      0, -6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), -ElasticModulus * CrossSecArea / (10 * a), 0, 6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), -ElasticModulus * CrossSecArea / (10 * a);
      0, -ElasticModulus * CrossSecArea / (10 * a), -2 * ElasticModulus * CrossSecArea / 15, 0, ElasticModulus * CrossSecArea / (10 * a), ElasticModulus * CrossSecArea / 30;
      0, 0, 0, 0, 0, 0;
      0, 6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), ElasticModulus * CrossSecArea / (10 * a), 0, -6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), ElasticModulus * CrossSecArea / (10 * a)
      0, -ElasticModulus * CrossSecArea / (10 * a), ElasticModulus * CrossSecArea / 30, 0, ElasticModulus * CrossSecArea / (10 * a), -2 * ElasticModulus * CrossSecArea / 15;
      ], [0, -3 * ElasticModulus * CrossSecArea / (5 * a ^ 2), -ElasticModulus * CrossSecArea / (10 * a), 0, 6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), -ElasticModulus * CrossSecArea / (10 * a);
      0, 0, 0, 6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), 0, 0;
      0, 0, 0, ElasticModulus * CrossSecArea / (10 * a), 0, 0;
      0, 3 * ElasticModulus * CrossSecArea / (5 * a ^ 2), ElasticModulus * CrossSecArea / (10 * a), 0, -6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), ElasticModulus * CrossSecArea / (10 * a);
      0, 0, 0, -6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), 0, 0;
      0, 0, 0, ElasticModulus * CrossSecArea / (10 * a), 0, 0;
      ], [0, 0, -ElasticModulus * CrossSecArea / 15, 0, ElasticModulus * CrossSecArea / (10 * a), ElasticModulus * CrossSecArea / 30;
      0, 0, 0, ElasticModulus * CrossSecArea / (10 * a), 0, 0;
      0, 0, 0, 2 * ElasticModulus * CrossSecArea / 15, 0, 0;
      0, 0, ElasticModulus * CrossSecArea / 15, 0, -ElasticModulus * CrossSecArea / (10 * a), -ElasticModulus * CrossSecArea / 30;
      0, 0, 0, -ElasticModulus * CrossSecArea / (10 * a), 0, 0;
      0, 0, 0, -ElasticModulus * CrossSecArea / 30, 0, 0
      ], [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, -6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), ElasticModulus * CrossSecArea / (10 * a);
      0, 0, 0, 0, -ElasticModulus * CrossSecArea / (10 * a), -ElasticModulus * CrossSecArea / 30;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 6 * ElasticModulus * CrossSecArea / (5 * a ^ 2), -ElasticModulus * CrossSecArea / (10 * a);
      0, 0, 0, 0, -ElasticModulus * CrossSecArea / (10 * a), 2 * ElasticModulus * CrossSecArea / 15
      ], [0, 0, 0, 0, -3 * ElasticModulus * CrossSecArea / (5 * a ^ 2), ElasticModulus * CrossSecArea / (10 * a);
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / (5 * a ^ 2), -ElasticModulus * CrossSecArea / (10 * a);
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0
      ], [0, 0, 0, 0, 0, -ElasticModulus * CrossSecArea / 15;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, ElasticModulus * CrossSecArea / 15;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0]]);

% Define stiffness matrix for one element
K3E = [sparse(6, 42), [0, 0, 0, 0, 0, 0;
      0, 36 * ElasticModulus * CrossSecArea / (35 * a ^ 3), 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0, - 108 * ElasticModulus * CrossSecArea / (35 * a ^ 3), 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2);
      0, 9 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 9 * ElasticModulus * CrossSecArea / (70 * a), 0, - 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0;
      0, 0, 0, 0, 0, 0;
      0, - 36 * ElasticModulus * CrossSecArea / (35 * a ^ 3), - 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0, 108 * ElasticModulus * CrossSecArea / (35 * a ^ 3), - 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2);
      0, 9 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0, 0, - 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 9 * ElasticModulus * CrossSecArea / (70 * a); ], [0, 0, 0, 0, 0, 0;
      0, 0, 9 * ElasticModulus * CrossSecArea / (70 * a), 0, - 27 * ElasticModulus * CrossSecArea / (35 * a ^ 2), 0;
      0, 0, - 3 * ElasticModulus * CrossSecArea / 280, 0, - 9 * ElasticModulus * CrossSecArea / (35 * a), 3 * ElasticModulus * CrossSecArea / 140;
      0, 0, 0, 0, 0, 0;
      0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * a), 0, 27 * ElasticModulus * CrossSecArea / (35 * a ^ 2), 0;
      0, 0, 3 * ElasticModulus * CrossSecArea / 280, 0, 0, 3 * ElasticModulus * CrossSecArea / 140; ], sparse(6, 6), [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 108 * ElasticModulus * CrossSecArea / (35 * a ^ 3), - 27 * ElasticModulus * CrossSecArea / (35 * a ^ 2);
      0, 0, 0, 0, 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, - 108 * ElasticModulus * CrossSecArea / (35 * a ^ 3), 27 * ElasticModulus * CrossSecArea / (35 * a ^ 2);
      0, 0, 0, 0, 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), - 9 * ElasticModulus * CrossSecArea / (35 * a); ], [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 9 * ElasticModulus * CrossSecArea / (70 * a);
      0, 0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / 280;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * a);
      0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 280; ], sparse(6, 12), [0, 0, 0, 0, 0, 0;
      0, 0, - ElasticModulus * CrossSecArea / 280, 0, - 9 * ElasticModulus * CrossSecArea / (70 * a), 3 * ElasticModulus * CrossSecArea / 280;
      0, 0, ElasticModulus * CrossSecArea * a / 35, 0, 3 * ElasticModulus * CrossSecArea / 280, - 3 * ElasticModulus * CrossSecArea * a / 280;
      0, 0, 0, 0, 0, 0;
      0, 0, ElasticModulus * CrossSecArea / 280, 0, 9 * ElasticModulus * CrossSecArea / (70 * a), - 3 * ElasticModulus * CrossSecArea / 280;
      0, 0, - ElasticModulus * CrossSecArea * a / 280, 0, - 3 * ElasticModulus * CrossSecArea / 280, ElasticModulus * CrossSecArea * a / 140; ], sparse(6, 6), [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0;
      0, 0, 0, 0, 9 * ElasticModulus * CrossSecArea / (70 * a), - 3 * ElasticModulus * CrossSecArea / 140;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, - 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0;
      0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 140; ], [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / 280;
      0, 0, 0, 0, 0, ElasticModulus * CrossSecArea * a / 140;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 280;
      0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea * a / 280; ], sparse(6, 60), [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, - 36 * ElasticModulus * CrossSecArea / (35 * a ^ 3), 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2);
      0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 0;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 36 * ElasticModulus * CrossSecArea / (35 * a ^ 3), - 27 * ElasticModulus * CrossSecArea / (70 * a ^ 2);
      0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * a ^ 2), 9 * ElasticModulus * CrossSecArea / (70 * a); ], [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * a);
      0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 280;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 9 * ElasticModulus * CrossSecArea / (70 * a);
      0, 0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / 280; ], sparse(6, 30), [0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, - ElasticModulus * CrossSecArea / 280;
      0, 0, 0, 0, 0, - ElasticModulus * CrossSecArea * a / 280;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, ElasticModulus * CrossSecArea / 280;
      0, 0, 0, 0, 0, ElasticModulus * CrossSecArea * a / 35; ]];


%% Assemble global system
% Updated this with the help of chatgpt to assemble more efficiently
% Preallocate cell arrays to store triplet data
Icell  = cell(nel, 1); Jcell  = cell(nel, 1);
VcellM = cell(nel, 1); VcellK = cell(nel, 1);
nvg=TotalDOFs;       % number of variables in global mesh

for ie = 1:nel
    % Extract global node numbers of the nodes of element ie
    % and the global coordinates of the element nodes
    nodes = mconn(ie,:);
    dofs = vec(DOFsPerNode*(nodes-1)+(1:DOFsPerNode).').';

    % Build triplets from local matrix
    [row_idx, col_idx] = ndgrid(dofs, dofs);
    
    Icell{ie} = row_idx(:); Jcell{ie} = col_idx(:);
    VcellM{ie} = M1E(:); VcellK{ie} = K1E(:);
    
end

% Concatenate all triplet data
I = vertcat(Icell{:}); J = vertcat(Jcell{:});

% Assemble linear global matrices
M1G  = sparse(I, J, vertcat(VcellM{:}), nvg, nvg);
K1G = sparse(I, J, vertcat(VcellK{:}), nvg, nvg);

% Assemble quadratic stiffness matrix 
Icell = cell(nel, 1); Jcell = cell(nel, 1); Vcell = cell(nel, 1);
for ie = 1:nel
    % Local node indices
    nodes = mconn(ie,:);
    dofs = vec(DOFsPerNode*(nodes-1)+(1:DOFsPerNode).').';
    dofsm1 = dofs - 1;
    dofs2 = vec(dofs' + dofsm1*nvg).';
    
    % Build triplets from local matrix
    [row_idx, col_idx] = ndgrid(dofs, dofs2);
    
    Icell{ie} = row_idx(:); Jcell{ie} = col_idx(:); Vcell{ie} = K2E(:);
end

% Concatenate all triplet data
I = vertcat(Icell{:}); J = vertcat(Jcell{:}); V = vertcat(Vcell{:});

% Assemble global sparse matrix
K2G = sparseCSR(I, J, full(V), nvg, nvg^2);

% Assemble cubic stiffness matrix 
Icell = cell(nel, 1); Jcell = cell(nel, 1); Vcell = cell(nel, 1);
for ie = 1:nel
    % Local node indices
    nodes = mconn(ie,:);
    dofs = vec(DOFsPerNode*(nodes-1)+(1:DOFsPerNode).').';
    dofsm1 = dofs - 1;
    dofs2 = vec(dofs' + dofsm1*nvg).';
    dofs3 = vec(dofs2' + dofsm1*nvg^2).';
    
    % Build triplets from local matrix
    [row_idx, col_idx] = ndgrid(dofs, dofs3);
    
    Icell{ie} = row_idx(:); Jcell{ie} = col_idx(:); Vcell{ie} = K3E(:);
end

% Concatenate all triplet data
I = vertcat(Icell{:}); J = vertcat(Jcell{:}); V = vertcat(Vcell{:});

% Assemble global sparse matrix
K3G = sparseCSR(I, J, full(V), nvg, nvg^3);


%% Impose boundary conditions
fixedDOFs = [1, 2, 3];
freeDOFs = setdiff(1:TotalDOFs, fixedDOFs);

% Reduced system method
K1G = K1G(freeDOFs, freeDOFs);
M1G = M1G(freeDOFs, freeDOFs);

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

D1G = d1 * M1G + d2 * K1G; % Add some damping for numerical stability
%% Convert to state-space representation
% At present, we have a second-order mechanical system with governing equations
%
%     M1G q̈ + D1G q̇ + K1G q + K2G (q ⊗ q) + K3G (q ⊗ q ⊗ q) + ...
%       = RB0 u
%     y = C₁ q̇ + C₂ q.
%
%   Defining the state vector x = [q q̇]^T, we can convert to a first-order
%   nonlinear state-space
%
%     ẋ = A x + N₂ (x ⊗ x) + N₃ (x ⊗ x ⊗ x) ...
%          + G₀ u + G₁ (x ⊗ u) + G₂ (x ⊗ x ⊗ u) + G₃ (x ⊗ x ⊗ x ⊗ u)
%     y  = C₁ q̇ + C₂ q.

      n = length(M1G);
      
      McholL = chol(M1G).'; % Use cholesky factor for inverting rather than inv()
      A = [sparse(n, n), speye(n);
            -McholL.' \ (McholL \ K1G), -McholL.' \ (McholL \ D1G)];
            
      % Construct F₂
      p = 2;
      idxs = vec(vec((1:n).' + (0:2 * n:2 * n * (n - 1))) + [0, (2 * n) ^ p / 2 + (2 * n) ^ (p - 1) / 2 + n * (p - 2)]);
      In2 = sparse(2 * n ^ p, (2 * n) ^ p);
      In2(:, idxs) = speye(2 * n ^ p);
      
      F2 = [sparse(n, n ^ 2), sparse(n, n ^ 2);
            -McholL.' \ (McholL \ K2G), sparse(n, n ^ 2)] * In2;
      
      % Construct F₃
      p = 3;
      idxs = vec(vec(vec((0:(2 * n) ^ (1 - 1):(2 * n) ^ (1 - 1) * (n - 1)).' + 1 + (0:(2 * n) ^ (2 - 1):(2 * n) ^ (2 - 1) * (n - 1))) + (0:(2 * n) ^ (3 - 1):(2 * n) ^ (3 - 1) * (n - 1))) + [0, (2 * n) ^ p / 2 + (2 * n) ^ (p - 1) / 2 + n * (p - 2)]);
      In3 = sparse(2 * n ^ p, (2 * n) ^ p);
      In3(:, idxs) = speye(2 * n ^ p);
      
      F3 = [sparse(n, n ^ 3), sparse(n, n ^ 3);
            -McholL.' \ (McholL \ K3G), sparse(n, n ^ 3)] * In3;
      
%% Format outputs
f = {A, F2, F3};
g = {eye(2*n)};
h = {eye(2*n)};

end

function [mconn,x] = generateMesh(a,nxe)
% Generate a uniform linear mesh for a linear domain
%
% Inputs: a - length of the domain
%       nxe - number of elements in x-direction
%
% Outputs: mconn - connectivity matrix
%            x - x locations for all of the nodes
%

nx=nxe+1;                               % number of nodes in each direction
x = linspace(0,a,nx);

% Construct connectivity matrix for the mesh (maps global nodes to element nodes)
mconn = zeros(nxe,2); % Initialize connectivity matrix
for i=1:nxe
    mconn(:,1)=1:nxe;
    mconn(:,2)=(1:nxe)+1;
end

end