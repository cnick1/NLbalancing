function [E, f, g, h, xg] = getSystem30(numElements, p, q)
%getSystem30  Generates a 1D nonlinear heat equation FEM model.
%   This model is an unsteady nonlinear heat equation (reaction-diffusion) on
%   a 1D domain with linear elements. The model is subject to Neumann boundary
%   control via the secondary variables, so there are 2 inputs. You can use just
%   one column of the B matrix to alternatively consider a model with 1
%   insulated boundaries and just one controlled boundary. The output is the
%   mean temperature. Goal is to bring the system to zero from some initial
%   condition.
%
%   Usage:   [f,g,h] = getSystem30()
%
%   Inputs: numElements    - number of elements to discretize each direction with
%                            (default = 19, leading to n=20 dimension model)
%           p - cell array of coefficients for the p(x,u) uₓₓ(x,t) term;
%               e.g. the diffusion coefficient
%           q - cell array of coefficients for the q(x,u) u(x,t) term;
%               e.g. the reaction coefficients
%           (p and q default to the Allen-Cahn equation if not supplied)
%
%   Outputs:    E - mass matrix for dynamics in generalized form
%          f,g,h  - Cell arrays containing the polynomial coefficients for
%                   the drift, input, and output in generalized form
%             xg - global node locations
%
%   Description: The underlying model is most similar to the PDE for
%   getSystem9(), but that model uses Chebychev discretization whereas here we
%   use finite elements. The model in getSystem27() is also similar: it is a 1D
%   unsteady nonlinear heat equation discretized with finite elements, but the
%   details of the PDE are slightly different so it is not an Allen-Cahn
%   equation. The governing PDE is
%
%     γ(x,u) uₜ(x,t) - p(x,u) uₓₓ(x,t) + q(x,u) u(x,t) = f(x,t)
%
%   subject to the boundary conditions
%       uₓ(0,t) = u₁(t)  (Neumann control input on end 1)
%       uₓ(1,t) = u₂(t)  (Neumann control input on end 2)
%
%   Any of the boundary conditions can be set to zero for an insulated boundary.
%   We assume γ(x,u)=1 for now and polynomial forms for p(x,u) and q(x,u):
%       p(x,u) = p₁(x) u + p₂(x) u^2 + p₃ u^3,
%       q(x,u) = q₁(x) u + q₂(x) u^2 + q₃ u^3,
%   so the function inputs are the cell arrays p = {p₁, p₂, p₃} and q = {q₁, q₂, q₃}.
%   After finite element discretization, the finite element equations for the
%   reaction-diffusion problem can be written as
%
%     M ẋ + K₁ x + K₃ (x⊗x⊗x) = B u,
%     y = C x.
%
%   which can of course be put in the standard state-space form by multiplying
%   by the inverse of the mass matrix. Here, we leave the mass matrix E=M in
%   place to preserve the sparsity of the model, which permits scaling the model
%   to fairly high dimensions.
%
%   For different p(x,u) and q(x,u), you can get different models. For example,
%   with p={ε,0,0} gives a standard linear diffusion operation. With q{1} = -λ,
%   you can model heat generation due to electrical resistance. With q{3} = -μ,
%   you may have a cubic source/reaction term such as in Allen-Cahn.
% TODO: Summarize function calls for different model types
%
%   The following table summarizes some options for the functions p(x,u) and
%   q(x,u) for a few different mathematical models. NEED TO CHECK SIGNS HERE
%
%                   Material parameters for different models
%   +-----------------------+-----------+-----------+-------------------------+
%   |      model name       |     p     |     q     |    brief description    |
%   +-----------------------+-----------+-----------+-------------------------+
%   |       heat eq.        |  {ε,0,0}  |  {0,0,0}  |  linear heat equation   |
%   |   unsteady heat eq.   |  {ε,0,0}  |  {λ,0,0}  |  wire w/ elec. resist.  |
%   |     Allen-Cahn eq.    |  {ε,0,0}  |  {1,0,-1} |    bistable equation    |
%   +-----------------------+-----------+-----------+-------------------------+
%
%   For example:
%       linear heat equation         - getSystem30(10, {0.1,0,0},  {0,0,0})
%       wire w/ elec resist from [3] - getSystem30(10, {1,0,0},    {3,0,0})
%       Allen-Cahn similar to [4]    - getSystem30(10, {0.01,0,0}, {1,0,-1})
%
%
%   Reference: [1] J. N. Reddy, An introduction to nonlinear finite element
%              analysis. Oxford University Press, 2004, doi:
%              10.1093/acprof:oso/9780198525295.001.0001
%              [2] S. A. Ragab and H. E. Fayed, Introduction to finite element
%              analysis for engineers. Taylor & Francis Group, 2017
%              [3] D. M. Boskovic, M. Krstic, and W. Liu, "Boundary control
%              of an unstable heat equation via measurement of
%              domain-averaged temperature,” IEEE Transactions on Automatic
%              Control, vol. 46, no. 12, pp. 2022–2028, 2001
%              [4] N. Corbin TBD
%
%   Part of the PPR repository.
%%

if nargin < 3
    q = {1,0,-1};
    if nargin < 2
        p = {0.01,0,0};
        if nargin < 1
            numElements = 19;
        end
    end
end

%% Prepare linear domain mesh parameters
nel = numElements;         % number of elements
nng = nel + 1;             % number of nodes global
nvpn = 1;                  % number of variables per node
nnpe = 2;                  % number of nodes per element
nvpe = nnpe * nvpn;        % number of variables per element
nvg = nng * nvpn;          % number of variables in global mesh
L = 1;                     % domain length

%% Specify material properties and forcing
% Material properties gamma, p(x,u), q(x,u)
gamma = 1; % mass property that has been nondimensionalized
% p and q are function inputs

% Construct forcing function f(x,y) (element-wise constant)
fo = 0; % Currently no force; leaving in case I want to modify it
fgen=zeros(nel,1);
igen=nel;
for ie=1:igen
    fgen(ie)=fo;
end

%% Generate a rectangular mesh on a rectangular domain
% Get connectivity matrix and node locations
[mconn,xg] = generateLinearMesh(L,nel);

%% Generate element matrices Me, Ke, Fe, and Re and Assemble global system Mg, Kg, Fg, Rg
[Mg,K1g,K2g,K3g,Fg,Rg] = assembleGlobalSystem(gamma,p,q,fgen,mconn,xg);

%% Prescribe and apply boundary conditions
% Here we will impose Neumann boundary conditions, which simply prescribe the secondary
% variable on the RHS, so the mass and stiffness are not complicated. Dirichlet are also
% possible, but they are a bit more complicated.
[Rg] = applyBoundaryConditions(Rg);

%% Put into standard control system form
% The current ODE we have for the global system is
%       Mg u̇ + K1g u + K2g (u ⊗ u) + K3g (u ⊗ u ⊗ u) = Fg + Rg.
% For computational purposes, we wish to put this in the descriptor state-space
% form for polynomial control systems:
%       Eẋ = A x + F₂ (x ⊗ x) + F₃ (x ⊗ x ⊗ x) + B u
%           y = C x
% made these negative at element level, saves a lot of time
A = K1g; % = -K1g; ^
F2 = K2g;
F3 = K3g;

% g(x) is going to be a B matrix; it is going to be computed by taking the F/R
% matrix and breaking up into two components, one for each boundary. Note,
% currently Fg is zero; the control input u is really going through the qs
% variables in Rg, so if you add "forcing" fgen later this would have to change.

B = sparse(nvg,2); % m=2 control inputs, one per boundary
B(1,1) = Rg(1);
B(end,2) = Rg(end);

%% Format outputs
E = Mg;
f = {A, F2, F3};
g = {B};
h = {1/nvg * ones(1,nvg)}; % average over domain


end

%% Helper functions

function [Rg] = applyBoundaryConditions(Rg)
% Prescribe and subsequently apply boundary conditions; here they are
% hard-coded (extra stuff was vestigial from the more complicated 2D example)

%% Prescribe boundary conditions using us and qs variables
% These will end up forming the Rg vector, which will be broken into B*u
% for 4 boundary control inputs.
% Dirichlet BCs specify u=us (none here)
%% Neumann BCs specify secondary variable
% hard-coded

% Prescribe boundary condition on end 1 (Neumann)
Rg(1) = 1; % # is just a placeholder, will be controlled by u

% Prescribe boundary condition on end 2 (Neumann)
Rg(end) = 1; % # is just a placeholder, will be controlled by u

end

function [Mg,K1g,K2g,K3g,Fg,Rg] = assembleGlobalSystem(gamma,p,q,fgen,mconn,xg)
% Generate element matrices and assemble the global system
% Note: At this stage no BCs will be applied yet
vec = @(X) X(:);
[nel,nnpe] = size(mconn); % number of elements & number of nodes per element
nng=length(xg);          % number of global nodes
nvpn=1;             % number of variables per node
nvpe=nnpe*nvpn;     % number of variables per element
nvg=nng*nvpn;       % number of variables in global mesh

%% Assemble cubic global matrix
% Compute element matrices outside of loop = assume same for every element
nodes = mconn(1,:); % get first element
xe = xg(nodes,1);
[Me,K1e,K2e,K3e,Fe] = mekefe(gamma,p,q,fgen(1),xe);
Re=zeros(nvpe,1);

%% Assemble linear global matrix
% Preallocate cell arrays to store triplet data
Icell  = cell(nel, 1); Jcell  = cell(nel, 1);
VcellM = cell(nel, 1); VcellK = cell(nel, 1);

% Preallocate global force and reaction vectors
Fg = zeros(nvg,1); Rg = zeros(nvg,1);

for ie = 1:nel
    % Extract global node numbers of the nodes of element ie
    % and the global x and y-coordinates of the element nodes
    nodes = mconn(ie,:);
    % xe = xg(nodes,1); ye = xg(nodes,2);
    
    %% Compute element mass Me, stiffness Ke, forcing Fe, and Re
    % If material property (P and q) and/or force f are element dependent,
    % we need to specify these properties here.
    % [Me,K1e,K3e,Fe] = mekefe(gamma,P,q,q{3},fgen(1),xe,ye);
    % Re=zeros(nvpe,1);
    
    % Build triplets from local matrix
    [row_idx, col_idx] = ndgrid(nodes, nodes);
    
    Icell{ie} = row_idx(:); Jcell{ie} = col_idx(:);
    VcellM{ie} = Me(:); VcellK{ie} = -K1e(:);
    
    Fg(nodes) = Fg(nodes) + Fe;
    Rg(nodes) = Rg(nodes) + Re;
end

% Concatenate all triplet data
I = vertcat(Icell{:}); J = vertcat(Jcell{:});

% Assemble global sparse matrix
% Since these are not short/wide, no need for sparseCSR or sparseIJV
Mg  = sparse(I, J, vertcat(VcellM{:}), nvg, nvg);
K1g = sparse(I, J, vertcat(VcellK{:}), nvg, nvg);

%% Assemble quadratic global matrix
% Preallocate cell arrays to store triplet data
Icell = cell(nel, 1);
Jcell = cell(nel, 1);
Vcell = cell(nel, 1);

for ie = 1:nel
    % assume K2e is constant (it was already computed)
    
    % Local node indices
    nodes = mconn(ie, :);
    nodesm1 = nodes - 1;
    nodes2 = vec(nodes' + nodesm1*nvg).';
    
    % Build triplets from local matrix
    [row_idx, col_idx] = ndgrid(nodes, nodes2);
    
    Icell{ie} = row_idx(:);
    Jcell{ie} = col_idx(:);
    Vcell{ie} = -K2e(:);
end

% Concatenate all triplet data
I = vertcat(Icell{:});
J = vertcat(Jcell{:});
V = vertcat(Vcell{:});

% Assemble global sparse matrix
K2g = sparseCSR(I, J, full(V), nvg, nvg^2);

%% Assemble cubic global matrix
% Preallocate cell arrays to store triplet data
Icell = cell(nel, 1);
Jcell = cell(nel, 1);
Vcell = cell(nel, 1);

for ie = 1:nel
    % [~,~,K3e,~] = mekefe(gamma,P,q,q{3},fgen(ie),xe,ye);
    % commented out = assume K3e is constant (it was already computed)
    
    % Local node indices
    nodes = mconn(ie, :);
    nodesm1 = nodes - 1;
    nodes2 = vec(nodes' + nodesm1*nvg).';
    nodes3 = vec(nodes2' + nodesm1*nvg^2).';
    
    % Build triplets from local matrix
    [row_idx, col_idx] = ndgrid(nodes, nodes3);
    
    Icell{ie} = row_idx(:);
    Jcell{ie} = col_idx(:);
    Vcell{ie} = -K3e(:);
end

% Concatenate all triplet data
I = vertcat(Icell{:});
J = vertcat(Jcell{:});
V = vertcat(Vcell{:});

% Assemble global sparse matrix
% Either sparseCSR or sparseIJV can be used (sparseCSR is actually more
% memory efficient). For large nvg though, Matlab won't permit using
% sparseCSR (which uses the transpose of sparse), so sparseIJV is necessary
% even though it is actually slightly worse
% K3g = sparseIJV(I, J, full(V), nvg, nvg^3);
K3g = sparseCSR(I, J, full(V), nvg, nvg^3);

end

function [mconn,xg] = generateLinearMesh(L,nel)
% Generate a uniform linear mesh for a linear domain
%
% Inputs: L - length of the domain
%       nel - number of elements
%
% Outputs: mconn - connectivity matrix
%             xg - x locations for all of the nodes
%

nng=nel+1;                 % total number of nodes
xg = linspace(0, L, nng).';  % node locations

% Construct connectivity matrix for the mesh (maps global nodes to element nodes)
mconn = zeros(nel,2); % Initialize connectivity matrix
mconn(:,1) = 1:nel; mconn(:,2) = 2:nel+1;

end

function [Me,K1e,K2e,K3e,Fe] = mekefe(gamma,p,q,f,xe)
% Return stiffness matrix and force vector (bilinear rectangular element)
% This function assumes elementwise constant coefficients p, q, and f.
%
% Inputs: gamma - mass property
%             p - material property
%             q - heat convection
%             f - distributed heat source
%            xe - x coordinates of element nodes; just for computing
%                 element widths
%
% Outputs:  Me  - element mass matrix
%           K1e - element linear stiffness matrix
%           K2e - element quadratic stiffness matrix
%           K3e - element cubic stiffness matrix
%           Fe  - element force vector
%

he=xe(2)-xe(1);

% Define element matrix components
L0 = he/6 * [2, 1;
    1, 2];

L1 = 1/he * [1, -1;
    -1,  1];

L2m = 1/2/he * [1, 0, 0, -1;
    -1, 0, 0, 1];
L2k = he/12 * [3, 2, 0, 1;
    1, 0, 2, 3];

L3m = 1/3/he * [1, 0, 0, 0, 0, 0, 0, -1;
    -1, 0, 0, 0, 0, 0, 0, 1];
L3k = he/20 * [4, 3, 0, 2, 0, 0, 0, 1;
    1, 0, 0, 0, 2, 0, 3, 4];

% Define element mass and stiffness matrices
Me = gamma*L0;
K1e = p{1} * L1  + q{1} * L0;
K2e = p{2} * L2m + q{2} * L2k;
K3e = p{3} * L3m + q{3} * L3k;

Fe = f*he/2*[1; 1];
end




