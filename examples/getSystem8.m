function [f, g, h] = getSystem8(numElements, lambda)
%getSystem8  Generates a cubic finite element heat equation model system
%   The system is a finite element model for a nonlinear heat equation,
%   i.e. a reaction-diffusion equation. The function returns a finite
%   element model with linear elements. (Note: A previous version used to
%   have quadratic elements partially implemented.) The state-space
%   dimension is n=numElements-1 due to the fixed end boundary conditions.
%
%   Usage:   [f,g,h] = getSystem8()
%
%   Inputs:
%       numElements    - number of elements to discretize the domain with
%                                                           (default =  4 )
%       lambda         - coefficient on the reaction term; drives
%                        instability for lambda > 0         (default = 1/8)
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The reaction-diffusion model used in [1] from [2-4] is a
%   represented by the nonlinear heat equation PDE
%
%     uₜ(x,t) = uₓₓ(x,t) + uₓ(x,t) + λ u(x,t) + u(x,t)³
%     u(0,t) = 0                     (Dirichlet BC)
%     u(1,t) = 0                     (Dirichlet BC)
%
%   which we augment with a control input. After FEM discretization, the
%   finite element equations can be written as
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
%   Reference: [1] N. A. Corbin and B. Kramer, “Scalable computation of
%               𝓗∞ energy functions for polynomial drift nonlinear
%               systems,” in 2024 American Control Conference (ACC), Jul.
%               2024, pp. 2506–2511. doi: 10.23919/acc60939.2024.10644363
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
%   Part of the PPR repository.
%%
vec = @(X) X(:);
if nargin < 2
    lambda = 1/8;
    if nargin < 1
        numElements = 4;
    end
end

%% Prepare linear domain mesh parameters
nel = numElements;            % number of elements
nng = nel + 1;                % number of nodes
nvpn = 1;                     % number of variables per node
nnpe = 2;                     % number of nodes per element
nvpe = nnpe * nvpn;           % number of variables per element
nvg = nng * nvpn;             % number of variables in global mesh

%% Specify material and geometry properties
% Define geometry and material properties
% I am using notation for the model problem (5.1) in Ragab [2]
L = 30;                       % wire length
gamma = 1;                    % mass property
p = 1;                        % nondimensionalized diffusion coefficient
q = -lambda;                  % -λ, to do with heat generation due to electrical resistance

%% Generate linear "mesh"
xg = linspace(0, L, nng);     % node locations
he = xg(nnpe) - xg(1);        % element length

%% Generate element matrices Me, Ke, and Re and Assemble global system Mg, Kg, Rg
% Define element matrix components
L0 = he/6 * [2, 1;
    1, 2];

L1 = 1/he * [1, -1;
    -1, 1];

L2 = 1/2*[-1, 1;
    -1, 1]; % If i remember correctly this is for the uₓ term

% Define element mass and stiffness matrices
Me = gamma*L0;
K1e = p*L1 + q*L0 + L2;

% Initialize and stack/assemble global matrix
Mg = zeros(nvg, nvg);
K1g = zeros(nvg, nvg);
for ie = 1:nel
    nodes = ie:(ie + nvpe - 1);
    
    % Assemble element matrices Me, Ke into global matrix Mg, Kg
    Mg(nodes, nodes)  = Mg(nodes, nodes)  + Me;
    K1g(nodes, nodes) = K1g(nodes, nodes) + K1e;
end


%% Assemble quadratic global matrix
% Form all zero quadratic global matrix directly
K2g = sparse(nvg, nvg^2);

%% Assemble cubic global matrix
% Define stiffness matrix for one element
K3e = -he/20 * [4, 3, 0, 2, 0, 0, 0, 1;
    1, 2, 0, 3, 0, 0, 0, 4];

% Initialize and stack/assemble global matrix
K3g = sparse(nvg, nvg^3);

for ie = 1:nel
    i = ie - 1;
    % Construct nodes and cnodes: nodes contains the indices for the
    % element nodes w.r.t. the global nodes. This is to map the element
    % indices to the global variable u. Similarly, nodes3 maps the indices
    % for the element nodes to the cubic variable (u⊗u⊗u). nodes3 is in
    % particular a little tricky, as we need to figure out:
    %   1) the skip due to going from element 1, 2, 3, i.e. the skip
    %      corresponding to ie increasing in the loop
    %       (for the linear variable, we just go up by 1 each time,
    %        but for the cubic variable we need to go up by chunks)
    %   2) the skips due to global nodes that are not in the element
    %              (for the linear variable, all the
    %              element nodes are grouped together)
    % Both of these skips are somewhat recursive, in that the skips for
    % (u⊗u⊗u) involve the skips needed for (u⊗u), etc. A commonly used
    % trick for this process is to add a column vector and a row vector to
    % get all the combinations.
    nodes = ie:(ie + nvpe - 1); % node numbers in u
    
    nodes3 = (nvg^2 + nvg + 1) * i  ... % shift (1) in starting index for on element iteration
        + vec(vec([1 2]' + [0 nvg]) ... % shift (2); [0 nvg] part basically does the quadratic
        + [0 nvg^2]);                   % shifts, which is then added to [0 nvg^2] for the cubic shifts
    
    % "stack" element matrices into global matrix
    K3g(nodes, nodes3) = K3g(nodes, nodes3) + K3e;
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

RB0 = generate_B_matrix(4, (nvg - 1) / 4);

%% Impose boundary conditions
fixedDOFs = [1, nvg]; % First and last nodes fixed
freeDOFs = setdiff(1:nvg, fixedDOFs);

% Reduced system method
K1g = K1g(freeDOFs, freeDOFs);
Mg  =  Mg(freeDOFs, freeDOFs);
RB0 = RB0(freeDOFs, :);

% K2G could clean up
fixedDOFsSquared = vec(fixedDOFs.' + (0:nvg:nvg ^ 2 - 1));
fixedDOFsSquared = unique([fixedDOFsSquared; vec((1:nvg).' + (fixedDOFs - 1) * nvg)]);

freeDOFsSquared = setdiff(1:nvg ^ 2, fixedDOFsSquared);
K2g = K2g(freeDOFs, freeDOFsSquared);

% K3G could clean up
fixedDOFsCubed = vec(vec((fixedDOFs - 1) * nvg + [1:nvg].') + (0:nvg ^ 2:nvg ^ 3 - 1));
fixedDOFsCubed = [fixedDOFsCubed; vec((1:nvg ^ 2).' + (fixedDOFs - 1) * nvg ^ 2)]; % Top rows zero
fixedDOFsCubed = [fixedDOFsCubed; vec(vec(fixedDOFs.' + (0:nvg:nvg ^ 2 - 1)) + (0:nvg ^ 2:nvg ^ 3 - 1))];
fixedDOFsCubed = sort(unique(fixedDOFsCubed));

freeDOFsCubed = setdiff(1:nvg ^ 3, fixedDOFsCubed);
K3g = K3g(freeDOFs, freeDOFsCubed);

%% Convert to state-space representation
n = length(Mg);

Mchol = chol(Mg).'; % Use Cholesky factor for inverting rather than inv()

% Construct F₁, F₂, & F₃
F1 = -Mchol.' \ (Mchol \ K1g);
F2 = -Mchol.' \ (Mchol \ K2g);
F3 = -Mchol.' \ (Mchol \ K3g);

% Construct B & C
B = Mchol.' \ (Mchol \ RB0);
C = RB0.';



%% Format outputs
f = {full(F1), F2, F3};
g = {full(B)};
h = {full(C)};

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
