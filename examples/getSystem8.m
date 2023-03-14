function [A, B, C, N, Q, f, g, h] = getSystem8(numElements, actuatorConfiguration, rotaryInertiaBool, varargin)
%getSystem8  Generates a cubic system for testing energy functions.
%            The system is a finite element model for a nonlinear (due to
%            von Karman strains) Euler-Bernoulli Beam.
%
%   Usage:  [A,B,C,N,Q] = getSystem8()
%        or [~,~,~,~,~,f,g,h] = getSystem8()
%%

% TODO: Code up true angle dependent bilinear inputs Q

vec = @(X) X(:);

if nargin < 1
  numElements = 3;
end

if nargin < 2
  actuatorConfiguration = 1;
end

if nargin < 3
  rotaryInertiaBool = false;
end



%% Define beam geometry and properties
BeamLength = 1; % length of beam
ElasticModulus = 210e9; % Young's modulus
CrossSecArea = 1e-1; % cross-sectional area
MomOfInertia = 1e-2; % moment of inertia
density = 8000; % density
delta = 0.1; % Cable attachment distance from centerline of beam if actuatorConfiguration = 2

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
  + rotaryInertiaBool*density * MomOfInertia / (30 * elementLength) * ...
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
switch actuatorConfiguration
    case 1
        RB0 = sparse(TotalDOFs, 2);        
        RB1 = sparse(TotalDOFs, 2*TotalDOFs);
        RB2 = sparse(TotalDOFs, 2*TotalDOFs^2); 
        RB3 = sparse(TotalDOFs, 2*TotalDOFs^3);  
        
        RB0(TotalDOFs - 2, 2) = 1; % Force in x direction
        RB0(TotalDOFs - 1, 1) = 1; % Force in y direction
        RB0(TotalDOFs, :) = 0;     % Moment in z direction
        
    case 2
        RB0 = sparse(TotalDOFs, 2);
        RB1 = sparse(TotalDOFs, 2*TotalDOFs);
        RB2 = sparse(TotalDOFs, 2*TotalDOFs^2); 
        RB3 = sparse(TotalDOFs, 2*TotalDOFs^3);
        
        % x     or u     -> (TotalDOFs - 2)
        % y     or v     -> (TotalDOFs - 1)
        % theta or dv/dx -> (TotalDOFs)
        
        % Order-0 terms
        RB0(TotalDOFs-2, 1) = -1;     % Force in x direction from cable 1
        RB0(TotalDOFs-2, 2) = -1;     % Force in x direction from cable 2
        RB0(TotalDOFs-1, 1) =  0;     % Force in y direction from cable 1
        RB0(TotalDOFs-1, 2) =  0;     % Force in y direction from cable 2
        RB0(TotalDOFs  , 1) =  delta; % Moment in z direction from cable 1
        RB0(TotalDOFs  , 2) = -delta; % Moment in z direction from cable 2
        
        % Order-1 terms (only appear in transverse direction)
        RB1(TotalDOFs-1, 2*(TotalDOFs-1)-1) = -1; % Force in y direction from cable 1 xn-1 u1
        RB1(TotalDOFs-1, 2*(TotalDOFs-1))   = -1; % Force in y direction from cable 2 xn-1 u2
        
        % Order-2 terms
        RB2(TotalDOFs-2, 2*(TotalDOFs-1)^2-1)             =  1/(2*L^2);         % Force in x direction from cable 1
        RB2(TotalDOFs-2, 2*(TotalDOFs-1)^2)               =  1/(2*L^2);         % Force in x direction from cable 2
        RB2(TotalDOFs-1, 2*(TotalDOFs-1)*(TotalDOFs-2)-1) =  1/L^2;             % Force in y direction from cable 1 
        RB2(TotalDOFs-1, 2*(TotalDOFs-1)*(TotalDOFs-2))   =  1/L^2;             % Force in y direction from cable 2 
        RB2(TotalDOFs  , 2*(TotalDOFs-1)^2-1)             = -delta/(2*L^2);     % Moment in z direction from cable 1 
        RB2(TotalDOFs  , 2*(TotalDOFs-1)^2)               =  delta/(2*L^2);     % Moment in z direction from cable 2 
        
        % x1x1u1 x1x1u2 x1x2u1 x1x2u2 x1x3u1 x1x3u2 ..... xnx1u1 xnx1u2 ... xnxnu1 xnxnu2
        % Need (xn-1)(xn-1) u1 = -1
        %      (xn-1)(xn-1) u2 = -1
        %
        % and (xn-1)(xn-2) u1 = -1
        %     (xn-1)(xn-2) u2 = -1
        %                               (or (xn-2)(xn-1) u1 = -1
        
        % Order-3 terms TODO !!! NOT CORRECT
        RB3(TotalDOFs-2, 2*(TotalDOFs-2)*(TotalDOFs-1)^2-1) = -1/L^3;         % Force in x direction from cable 1
        RB3(TotalDOFs-2, 2*(TotalDOFs-2)*(TotalDOFs-1)^2)   = -1/L^3;         % Force in x direction from cable 2
        RB3(TotalDOFs-1, 2*(TotalDOFs-1)*(TotalDOFs-2)^2-1) = -1/L^3;         % Force in y direction from cable 1 part1
        RB3(TotalDOFs-1, 2*(TotalDOFs-1)*(TotalDOFs-2)^2)   = -1/L^3;         % Force in y direction from cable 2 part1
        RB3(TotalDOFs-1, 2*(TotalDOFs-1)^3-1)               =  1/(2*L^3);     % Force in y direction from cable 1 part2
        RB3(TotalDOFs-1, 2*(TotalDOFs-1)^3)                 =  1/(2*L^3);     % Force in y direction from cable 2 part2
        RB3(TotalDOFs  , 2*(TotalDOFs-2)*(TotalDOFs-1)^2-1) =  delta/L^3;     % Moment in z direction from cable 1 
        RB3(TotalDOFs  , 2*(TotalDOFs-2)*(TotalDOFs-1)^2)   = -delta/L^3;     % Moment in z direction from cable 2 
end
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
RB0 = RB0(freeDOFs, :);
RB1 = RB1(freeDOFs, [freeDOFs, freeDOFs+TotalDOFs]);

% K2G could clean up
fixedDOFsSquared = vec(fixedDOFs.' + (0:TotalDOFs:TotalDOFs^2-1));
fixedDOFsSquared = unique([fixedDOFsSquared; vec((1:TotalDOFs).' + (fixedDOFs-1)*TotalDOFs)]);

freeDOFsSquared = setdiff(1:TotalDOFs^2, fixedDOFsSquared);
K2G = K2G(freeDOFs,freeDOFsSquared);
RB2 = RB2(freeDOFs, [freeDOFsSquared, freeDOFsSquared+TotalDOFs^2]);



% K3G could clean up
fixedDOFsCubed = vec(vec((fixedDOFs-1)*TotalDOFs+[1:TotalDOFs].') + (0:TotalDOFs^2:TotalDOFs^3-1));
fixedDOFsCubed = [fixedDOFsCubed; vec((1:TotalDOFs^2).' + (fixedDOFs-1)*TotalDOFs^2)]; % Top rows zero
fixedDOFsCubed = [fixedDOFsCubed; vec(vec(fixedDOFs.' + (0:TotalDOFs:TotalDOFs^2-1)) + (0:TotalDOFs^2:TotalDOFs^3-1))];
fixedDOFsCubed = sort(unique(fixedDOFsCubed));

freeDOFsCubed = setdiff(1:TotalDOFs^3, fixedDOFsCubed);
K3G = K3G(freeDOFs,freeDOFsCubed);
RB3 = RB3(freeDOFs, [freeDOFsCubed, freeDOFsCubed+TotalDOFs^3]);


D1G = 0.00001 * M1G + 0.00001 * K1G;

%% Convert to state-space representation
n = length(M1G);

Minv = inv(M1G); % Yikes!

McholL = chol(M1G).';
Minv * K1G-McholL.'\(McholL\K1G) 

N1 = [sparse(n, n), speye(n);
     -Minv * K1G, -Minv * D1G];

G0 = [sparse(n, 2);
     Minv * RB0];

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
     -Minv * K3G, sparse(n,n^3)]*In3;
 
% Construct Q
Im = speye(2);
G1 = [sparse(n,2*n), sparse(n,2*n);
     Minv * RB1, sparse(n,2*n)];
 
G2 = [sparse(n,2*n^2), sparse(n,2*n^2);
     Minv * RB2, sparse(n,2*n^2)]*kron(In2,Im);
 
G3 = [sparse(n,2*n^3), sparse(n,2*n^3);
     Minv * RB3, sparse(n,2*n^3)]*kron(In3,Im);
 


%% Format outputs
f = {N1,N2,N3}; 
g = {G0,G1,G2,G3};
h = {C};

A = full(N1);
B = full(G0);
C = full(C);
N = full(N2);
Q = G1;


end
