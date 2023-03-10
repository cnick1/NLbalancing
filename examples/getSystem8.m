function [A, B, C, N, Q, f, g, h] = getSystem8(numElements, varargin)
%getSystem8  Generates a cubic system for testing energy functions.
%            The system is a finite element model for a nonlinear (due to
%            von Karman strains) Euler-Bernoulli Beam.
%
%   Usage:  [A,B,C,N,Q] = getSystem8()
%        or [~,~,~,~,~,f,g,h] = getSystem8()
%%

% TODO: Code up true angle dependent bilinear inputs Q

vec = @(X) X(:);

if nargin == 0
  numElements = 3;
end

%% Define beam geometry and properties
BeamLength = 1; % length of beam
ElasticModulus = 210e9; % Young's modulus
CrossSecArea = 1e-1; % cross-sectional area
MomOfInertia = 1e-2; % moment of inertia
density = 8000; % density

% Define element properties
numNodes = numElements + 1; % number of nodes
x = linspace(0, BeamLength, numNodes); % node locations
elementLength = x(2) - x(1); % element length

% Define DOF counts
DOFsPerNode = 3;
DOFsPerElement = 2 * DOFsPerNode;
TotalDOFs = numNodes * DOFsPerNode;

%% Assemble linear global matrices (mass and stiffness)
% Define mass matrix for one element
M1E = density * CrossSecArea * elementLength / 420 * ...
[140, 0, 0, 70, 0, 0;
 0, 156, 22 * elementLength, 0, 54, -13 * elementLength;
 0, 22 * elementLength, 4 * elementLength ^ 2, 0, 13 * elementLength, -3 * elementLength ^ 2;
 70, 0, 0, 140, 0, 0;
 0, 54, 13 * elementLength, 0, 156, -22 * elementLength;
 0, -13 * elementLength, -3 * elementLength ^ 2, 0, -22 * elementLength, 4 * elementLength ^ 2] ...
  + density * MomOfInertia / (30 * elementLength) * ...
  [0, 0, 0, 0, 0, 0;
 0, 36, 33 * elementLength, 0, -36, 3 * elementLength;
 0, 3 * elementLength, 4 * elementLength ^ 2, 0, -3 * elementLength, -elementLength ^ 2;
 0, 0, 0, 0, 0, 0;
 0, -36, -3 * elementLength, 0, 36, -33 * elementLength;
 0, 3 * elementLength, -elementLength ^ 2, 0, -3 * elementLength, 4 * elementLength ^ 2];

% Define stiffness matrix for one element
K1E = [ElasticModulus * CrossSecArea / elementLength, 0, 0, -ElasticModulus * CrossSecArea / elementLength, 0, 0;
       0, 12 * ElasticModulus * MomOfInertia / elementLength ^ 3, 6 * ElasticModulus * MomOfInertia / elementLength ^ 2, 0, -12 * ElasticModulus * MomOfInertia / elementLength ^ 3, 6 * ElasticModulus * MomOfInertia / elementLength ^ 2;
       0, 6 * ElasticModulus * MomOfInertia / elementLength ^ 2, 4 * ElasticModulus * MomOfInertia / elementLength, 0, -6 * ElasticModulus * MomOfInertia / elementLength ^ 2, 2 * ElasticModulus * MomOfInertia / elementLength;
       -ElasticModulus * CrossSecArea / elementLength, 0, 0, ElasticModulus * CrossSecArea / elementLength, 0, 0;
       0, -12 * ElasticModulus * MomOfInertia / elementLength ^ 3, -6 * ElasticModulus * MomOfInertia / elementLength ^ 2, 0, 12 * ElasticModulus * MomOfInertia / elementLength ^ 3, -6 * ElasticModulus * MomOfInertia / elementLength ^ 2;
       0, 6 * ElasticModulus * MomOfInertia / elementLength ^ 2, 2 * ElasticModulus * MomOfInertia / elementLength, 0, -6 * ElasticModulus * MomOfInertia / elementLength ^ 2, 4 * ElasticModulus * MomOfInertia / elementLength];

% Initialize and stack/assemble global matrix
M1G = sparse(TotalDOFs, TotalDOFs);
K1G = sparse(TotalDOFs, TotalDOFs);

for i = 0:numElements - 1 % start from zero so you don't have to subtract 1 every time
  ii = i * DOFsPerNode + 1;
  M1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) = M1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) + M1E;
  K1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) = K1G(ii:(ii + DOFsPerElement - 1), ii:(ii + DOFsPerElement - 1)) + K1E;
end

%% Assemble quadratic global matrix
% Define stiffness matrix for one element
K2E = sparse([[0, 0, 0, 0, 0, 0;
               0, -6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), -ElasticModulus * CrossSecArea / (10 * elementLength), 0, 6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), -ElasticModulus * CrossSecArea / (10 * elementLength);
               0, -ElasticModulus * CrossSecArea / (10 * elementLength), -2 * ElasticModulus * CrossSecArea / 15, 0, ElasticModulus * CrossSecArea / (10 * elementLength), ElasticModulus * CrossSecArea / 30;
               0, 0, 0, 0, 0, 0;
               0, 6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), ElasticModulus * CrossSecArea / (10 * elementLength), 0, -6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), ElasticModulus * CrossSecArea / (10 * elementLength)
               0, -ElasticModulus * CrossSecArea / (10 * elementLength), ElasticModulus * CrossSecArea / 30, 0, ElasticModulus * CrossSecArea / (10 * elementLength), -2 * ElasticModulus * CrossSecArea / 15;
               ], [0, -3 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), -ElasticModulus * CrossSecArea / (10 * elementLength), 0, 6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), -ElasticModulus * CrossSecArea / (10 * elementLength);
               0, 0, 0, 6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), 0, 0;
               0, 0, 0, ElasticModulus * CrossSecArea / (10 * elementLength), 0, 0;
               0, 3 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), ElasticModulus * CrossSecArea / (10 * elementLength), 0, -6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), ElasticModulus * CrossSecArea / (10 * elementLength);
               0, 0, 0, -6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), 0, 0;
               0, 0, 0, ElasticModulus * CrossSecArea / (10 * elementLength), 0, 0;
               ], [0, 0, -ElasticModulus * CrossSecArea / 15, 0, ElasticModulus * CrossSecArea / (10 * elementLength), ElasticModulus * CrossSecArea / 30;
               0, 0, 0, ElasticModulus * CrossSecArea / (10 * elementLength), 0, 0;
               0, 0, 0, 2 * ElasticModulus * CrossSecArea / 15, 0, 0;
               0, 0, ElasticModulus * CrossSecArea / 15, 0, -ElasticModulus * CrossSecArea / (10 * elementLength), -ElasticModulus * CrossSecArea / 30;
               0, 0, 0, -ElasticModulus * CrossSecArea / (10 * elementLength), 0, 0;
               0, 0, 0, -ElasticModulus * CrossSecArea / 30, 0, 0
               ], [0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, -6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), ElasticModulus * CrossSecArea / (10 * elementLength);
               0, 0, 0, 0, -ElasticModulus * CrossSecArea / (10 * elementLength), -ElasticModulus * CrossSecArea / 30;
               0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 6 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), -ElasticModulus * CrossSecArea / (10 * elementLength);
               0, 0, 0, 0, -ElasticModulus * CrossSecArea / (10 * elementLength), 2 * ElasticModulus * CrossSecArea / 15
               ], [0, 0, 0, 0, -3 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), ElasticModulus * CrossSecArea / (10 * elementLength);
               0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / (5 * elementLength ^ 2), -ElasticModulus * CrossSecArea / (10 * elementLength);
               0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 0, 0
               ], [0, 0, 0, 0, 0, -ElasticModulus * CrossSecArea / 15;
               0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 0, ElasticModulus * CrossSecArea / 15;
               0, 0, 0, 0, 0, 0;
               0, 0, 0, 0, 0, 0]]);

% Initialize and stack/assemble global matrix
K2G = sparse(TotalDOFs, TotalDOFs ^ 2);

for i = 0:numElements - 1 % start from zero so you don't have to subtract 1 every time
  ii = i * DOFsPerNode + 1;
  idxs = (TotalDOFs + 1) * i * DOFsPerNode ... % Starting index shift depending on element iteration
    + vec(( ...
    [1:DOFsPerElement] ... % [1,2,3,4,5,6] (basically the linear indices)
    + [0:TotalDOFs:TotalDOFs * (DOFsPerElement - 1)]' ... % Add skips into sequence (add row to column and then vec)
  )')';

  % "stack" element matrices into global matrix
  K2G(ii:(ii + DOFsPerElement - 1), idxs) = K2G(ii:(ii + DOFsPerElement - 1), idxs) + K2E;
end

%% Assemble cubic global matrix
% Define stiffness matrix for one element
K3E = [sparse(6, 42), [0, 0, 0, 0, 0, 0;
                       0, 36 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0, - 108 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2);
                       0, 9 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 9 * ElasticModulus * CrossSecArea / (70 * elementLength), 0, - 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0;
                       0, 0, 0, 0, 0, 0;
                       0, - 36 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), - 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0, 108 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), - 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2);
                       0, 9 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0, 0, - 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 9 * ElasticModulus * CrossSecArea / (70 * elementLength); ], [0, 0, 0, 0, 0, 0;
                       0, 0, 9 * ElasticModulus * CrossSecArea / (70 * elementLength), 0, - 27 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 2), 0;
                       0, 0, - 3 * ElasticModulus * CrossSecArea / 280, 0, - 9 * ElasticModulus * CrossSecArea / (35 * elementLength), 3 * ElasticModulus * CrossSecArea / 140;
                       0, 0, 0, 0, 0, 0;
                       0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * elementLength), 0, 27 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 2), 0;
                       0, 0, 3 * ElasticModulus * CrossSecArea / 280, 0, 0, 3 * ElasticModulus * CrossSecArea / 140; ], sparse(6, 6), [0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 108 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), - 27 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 2);
                       0, 0, 0, 0, 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0;
                       0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, - 108 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), 27 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 2);
                       0, 0, 0, 0, 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), - 9 * ElasticModulus * CrossSecArea / (35 * elementLength); ], [0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, 9 * ElasticModulus * CrossSecArea / (70 * elementLength);
                       0, 0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / 280;
                       0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * elementLength);
                       0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 280; ], sparse(6, 12), [0, 0, 0, 0, 0, 0;
                       0, 0, - ElasticModulus * CrossSecArea / 280, 0, - 9 * ElasticModulus * CrossSecArea / (70 * elementLength), 3 * ElasticModulus * CrossSecArea / 280;
                       0, 0, ElasticModulus * CrossSecArea * elementLength / 35, 0, 3 * ElasticModulus * CrossSecArea / 280, - 3 * ElasticModulus * CrossSecArea * elementLength / 280;
                       0, 0, 0, 0, 0, 0;
                       0, 0, ElasticModulus * CrossSecArea / 280, 0, 9 * ElasticModulus * CrossSecArea / (70 * elementLength), - 3 * ElasticModulus * CrossSecArea / 280;
                       0, 0, - ElasticModulus * CrossSecArea * elementLength / 280, 0, - 3 * ElasticModulus * CrossSecArea / 280, ElasticModulus * CrossSecArea * elementLength / 140; ], sparse(6, 6), [0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0;
                       0, 0, 0, 0, 9 * ElasticModulus * CrossSecArea / (70 * elementLength), - 3 * ElasticModulus * CrossSecArea / 140;
                       0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, - 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0;
                       0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 140; ], [0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / 280;
                       0, 0, 0, 0, 0, ElasticModulus * CrossSecArea * elementLength / 140;
                       0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 280;
                       0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea * elementLength / 280; ], sparse(6, 60), [0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, - 36 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2);
                       0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 0;
                       0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 36 * ElasticModulus * CrossSecArea / (35 * elementLength ^ 3), - 27 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2);
                       0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * elementLength ^ 2), 9 * ElasticModulus * CrossSecArea / (70 * elementLength); ], [0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, - 9 * ElasticModulus * CrossSecArea / (70 * elementLength);
                       0, 0, 0, 0, 0, - 3 * ElasticModulus * CrossSecArea / 280;
                       0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, 9 * ElasticModulus * CrossSecArea / (70 * elementLength);
                       0, 0, 0, 0, 0, 3 * ElasticModulus * CrossSecArea / 280; ], sparse(6, 30), [0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, - ElasticModulus * CrossSecArea / 280;
                       0, 0, 0, 0, 0, - ElasticModulus * CrossSecArea * elementLength / 280;
                       0, 0, 0, 0, 0, 0;
                       0, 0, 0, 0, 0, ElasticModulus * CrossSecArea / 280;
                       0, 0, 0, 0, 0, ElasticModulus * CrossSecArea * elementLength / 35; ]];

% Initialize and stack/assemble global matrix
K3G = sparse(TotalDOFs, TotalDOFs ^ 3);

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

%% RHS
RB = sparse(TotalDOFs, 2);
RB(TotalDOFs - 1, 1) = 1;
RB(TotalDOFs - 2, 2) = 1;

%% Impose boundary conditions
fixedDOFs = [1, 2, 3];
freeDOFs = setdiff(1:TotalDOFs, fixedDOFs);

% In-place method
% K1G(fixedDOFs, :) = 0; K1G(:, fixedDOFs) = 0;
% M1G(fixedDOFs, :) = 0; M1G(:, fixedDOFs) = 0;
% M1G(fixedDOFs, fixedDOFs) = speye(3);
% RB(fixedDOFs, :) = 0;

% Reduced system method
K1G = K1G(freeDOFs, freeDOFs);
M1G = M1G(freeDOFs, freeDOFs);
RB = RB(freeDOFs, :);

% K2G could clean up
fixedDOFsSquared = vec(fixedDOFs.' + (0:TotalDOFs:TotalDOFs^2-1));
fixedDOFsSquared = unique([fixedDOFsSquared; vec((1:TotalDOFs).' + (fixedDOFs-1)*TotalDOFs)]);

freeDOFsSquared = setdiff(1:TotalDOFs^2, fixedDOFsSquared);
K2G = K2G(freeDOFs,freeDOFsSquared);

% K3G could clean up
fixedDOFsCubed = vec(vec((fixedDOFs-1)*TotalDOFs+[1:TotalDOFs].') + (0:TotalDOFs^2:TotalDOFs^3-1));
fixedDOFsCubed = [fixedDOFsCubed; vec((1:TotalDOFs^2).' + (fixedDOFs-1)*TotalDOFs^2)]; % Top rows zero
fixedDOFsCubed = [fixedDOFsCubed; vec(vec(fixedDOFs.' + (0:TotalDOFs:TotalDOFs^2-1)) + (0:TotalDOFs^2:TotalDOFs^3-1))];
fixedDOFsCubed = sort(unique(fixedDOFsCubed));

freeDOFsCubed = setdiff(1:TotalDOFs^3, fixedDOFsCubed);
K3G = K3G(freeDOFs,freeDOFsCubed);


D1G = 0.00001 * M1G + 0.00001 * K1G;

%% Convert to state-space representation
n = length(M1G);

Minv = inv(M1G); % Yikes!

A = [sparse(n, n), speye(n);
     -Minv * K1G, -Minv * D1G];

B = [sparse(n, 2);
     Minv * RB];

C = sparse(1, 2 * n); C(1, n - 1) = 1;

% Construct N_2
p=2;
idxs = vec(vec((1:n).'+(0:2*n:2*n*(n-1))) + [0, (2*n)^p/2 + (2*n)^(p-1)/2+n*(p-2)]);
In2 = sparse(2*n^p, (2*n)^p);
In2(:,idxs) = speye(2*n^p);

N2 = [sparse(n, n^2), sparse(n, n^2);
     -Minv * K2G, sparse(n,n^2)]*In2;

% Construct N_3
p=3;
idxs = vec(vec(vec((0:(2*n)^(1-1):(2*n)^(1-1)*(n-1)).'+1 + (0:(2*n)^(2-1):(2*n)^(2-1)*(n-1)))+(0:(2*n)^(3-1):(2*n)^(3-1)*(n-1))) + [0, (2*n)^p/2 + (2*n)^(p-1)/2+n*(p-2)]);
In3 = sparse(2*n^p, (2*n)^p);
In3(:,idxs) = speye(2*n^p);

N3 = [sparse(n, n^3), sparse(n,n^3);
     -Minv * K3G, sparse(n,n^3)]*In3
 
% Construct Q
Q = 0;


%% Format outputs
f = {A,N2,N3}; 
g = {B,Q};
h = {C};

A = full(A);
B = full(B);
C = full(C);
N = full(N2);


end
