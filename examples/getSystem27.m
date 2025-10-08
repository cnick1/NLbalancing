function [f, g, h, xg] = getSystem27(numElements, eps, lambda, mu)
%getSystem27 Generates a FEM model for a 1D unstable heat equation for testing controllers.
%   The PDE described in [1] models a thin wire with electricity running
%   through it. The resistance of the wire causes internal heat generation,
%   which is counteracted by the radiation of heat from the wire. This
%   function gives a nonlinear generalization of the model in [1].
%
%   Usage:   [f,g,h] = getSystem27()
%
%   Inputs:
%       numElements    - number of elements to discretize the domain with
%                                                            (default = 19)
%       eps            - diffusion coefficient; has a stabilizing/smoothing
%                        effect for eps > 0                   (default = 1)
%       lambda         - coefficient on the linear reaction term; drives
%                        instability for lambda > 0           (default = 1)
%       mu             - coefficient on the cubic reaction term; produces
%                        stable wells away from 0 for mu < 0  (default =-1)
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%                  xg  - global node locations
%
%   Description: The general form for the PDE considered by this function is
%
%     uₜ(x,t)  = ε uₓₓ(x,t) + λ u(x,t) + μ u(x,t)³
%     uₓ(0,t) = 0                     (Neumann BC)
%     uₓ(1,t) = u(t)                  (Neumann boundary control input)
%
%   For different values of the parameters, different systems can be
%   obtained. For example, the original inspiration for this model is from
%   [1], which corresponds to ε=1, λ=3, and μ=0. This models a wire as heat
%   conduction in a thin rod of small constant cross section. All of the
%   material properties are taken as constant EXCEPT the electrical resistivity,
%   which varies linearly with respect to the temperature of the wire.
%
%   Another example that can be obtained is the Allen-Cahn equation can be
%   obtained (with Neumann rather than the usual Dirichlet boundary
%   conditions) with ε=1/100, λ=1, and μ=-3. The general form can be
%   considered a reaction-diffusion heat equation model, with one insulated
%   end and a Neumann control at the other end.(Nuemann control is easy to
%   implement since it is just in the secondary variable). After finite
%   element discretization, the finite element equations for the reaction-
%   diffusion problem can be written as
%
%     M ẋ + K₁ x + K₃ (x⊗x⊗x) = B u,
%     y = C x.
%
%   which can of course be put in the standard state-space form by
%   multiplying by the inverse of the mass matrix. In [1], the main results
%   correspond to Dirichlet boundary control, so here we made some changes
%   to more easily implement it with Neumann control (mentioned in [1] but
%   not demonstrated).
%
%   Reference: [1] D. M. Boskovic, M. Krstic, and W. Liu, "Boundary control
%              of an unstable heat equation via measurement of
%              domain-averaged temperature,” IEEE Transactions on Automatic
%              Control, vol. 46, no. 12, pp. 2022–2028, 2001
%              [2] S. A. Ragab and H. E. Fayed, Introduction to finite
%              element analysis for engineers. Taylor & Francis Group, 2017
%
%   Part of the PPR repository.
%%
vec = @(X) X(:);
if nargin < 4
    mu = -1;
    if nargin < 3
        lambda = 1;
        if nargin < 2
            eps = 1;
            if nargin < 1
                numElements = 19;
            end
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

%% Specify material and geometry properties
% Define geometry and material properties
% I am using notation for the model problem (5.1) in Ragab [2]
L = 1;                     % wire length
gamma = 1;                 % mass property
p  = eps;                  %  ε, diffusion coefficient
q  = -lambda;              % -λ, to do with heat generation due to electrical resistance
q3 = -mu;                  % -μ, cubic source/reaction term

%% Generate linear "mesh"
xg = linspace(0, L, nng);  % node locations
he = xg(2) - xg(1);        % element length

%% Generate element matrices Me, K1e, and Re and Assemble global system Mg, Kg, Rg
% Define element matrix components
L0 = he/6 * [2, 1;
    1, 2];

L1 = 1/he * [1, -1;
    -1,  1];

% Define element mass and stiffness matrices
Me = gamma*L0;
K1e = p*L1 + q*L0;

% Initialize and stack/assemble global matrix
Mg = zeros(nvg, nvg);
K1g = zeros(nvg, nvg);
Fg = zeros(nvg, 1);
Rg = zeros(nvg, 1);
for ie = 1:nel
    nodes = ie:(ie + nvpe - 1);
    
    % Assemble element matrices Me, K1e into global matrix Mg, Kg
    Mg(nodes,nodes) = Mg(nodes, nodes) + Me;
    K1g(nodes,nodes) = K1g(nodes, nodes) + K1e;
end

%% Assemble quadratic global matrix
% Form all zero quadratic global matrix directly
K2g = sparse(nvg, nvg^2);

%% Assemble cubic global matrix
% Define stiffness matrix for one element
L3 = he/20 * [4, 3, 0, 2, 0, 0, 0, 1;
    1, 0, 0, 0, 2, 0, 3, 4];

K3e = q3 * L3;

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

%% Prescribe and apply boundary conditions
% x=0 end is insultated; corresponding node is #1
% Rg(1) = Rg(1) + 0; % Rg is already zero

% x=L end has Neumann boundary control
Rg(end) = Rg(end) + 1; % 1 is just a placeholder for B; signal will be u(t)

%% Put into standard control system form
% The current ODE we have for the global system is Mg u̇ + Kg u = Fg + Rg.
% We wish to put this in the standard state-space form for control systems:
%       ẋ = A x + B u
%       y = C x
% Multiplying my M⁻¹ achieves this. The new system will be
%       ̇u = -Mg⁻¹Kg u + Mg⁻¹Fg + Mg⁻¹Rg
%        := A u  + F
% For computational purposes, we will do this inversion by computing the
% Cholesky factor of M; if we were to do M\K, M\F, and M\R these would all
% solve using the Cholesky factor, so we are recycling it. (We could
% simulate faster using sparse matrices and specifying mass matrix as an
% ode option; the inversion destroys the sparsity.)

Mchol = chol(Mg).'; % Use Cholesky factor for inverting
A = -Mchol.' \ (Mchol \ K1g);
% F2 = -Mchol.' \ (Mchol \ K2g);
F2 = sparse(nvg, nvg^2);
F3 = sparse(-Mchol.' \ (Mchol \ K3g));

% Now, changing notation a bit, we want to put the system in the form
%       ̇x = f(x) + g(x) u
% so the old u becomes x and we introduce this new u for controls. f(x) is
% easy, it is just f(x) = A x. g(x) is going to be a B matrix; it is going
% to be computed by taking the F matrix; the control input u is really
% going through the qs variables in Rg, so if you add "forcing" fgen later
% this would have to change.
B = Rg;
B =  Mchol.' \ (Mchol \ B);



%% Format outputs
f = {A, F2, F3};
g = {B};
h = {1/nvg * ones(1,nvg)}; % average over domain


end

