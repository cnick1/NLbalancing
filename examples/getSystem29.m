function [E, f, g, h, xyg] = getSystem29(numElements, eps, lambda, mu)
%getSystem29  Generates a FEM 2D Allen-Cahn equation model. This model is an
%   unsteady nonlinear heat equation (reaction-diffusion) on a square 2D domain.
%   The function returns a finite element model with bilinear rectangular
%   elements. The model is subject to Neumann boundary control via the secondary
%   variables, so there are 4 inputs. You can use just one column of the B
%   matrix to alternatively consider a model with 3 insulated boundaries and
%   just one controlled boundary. The output is the mean temperature. Goal is to
%   bring the system to zero from some initial condition.
%
%   Usage:   [f,g,h] = getSystem29()
%
%   Inputs: numElements    - number of elements to discretize each direction with
%                            (default = 8, leading to n=81 dimension model)
%           eps            - diffusion coefficient; has a stabilizing/smoothing
%                            effect for eps > 0 (default = 1)
%           lambda         - coefficient on the linear reaction term; drives
%                            instability for lambda > 0 (default = 1)
%           mu             - coefficient on the cubic reaction term; produces stable
%                            wells away from 0 for mu < 0  (default =-1)
%
%   Outputs:    E - mass matrix for dynamics in generalized form
%          f,g,h  - Cell arrays containing the polynomial coefficients for
%                   the drift, input, and output in generalized form
%             xyg - global node locations
%
%   Description: The underlying model is most similar to the PDE for
%   getSystem9(), but that model uses Chebychev discretization whereas here we
%   use finite elements. The model in getSystem27() is also similar: it is a 1D
%   unsteady nonlinear heat equation discretized with finite elements, but the
%   details of the PDE are slightly different so it is not an Allen-Cahn
%   equation. The governing PDE is
%
%     uₜ(x,y,t) = ε uₓₓ(x,y,t) + u(x,y,t) - u(x,y,t)³
%
%   subject to the boundary conditions
%       uᵧ(x,0,t) = u₁(t)  (Neumann control input on side AB)
%       uₓ(1,y,t) = u₂(t)  (Neumann control input on side BC)
%       uᵧ(x,1,t) = u₃(t)  (Neumann control input on side CD)
%       uₓ(0,y,t) = u₄(t)  (Neumann control input on side DA)
%
%   Any of the boundary conditions can be set to zero for an insulated boundary.
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
%   Reference: [1] J. N. Reddy, An introduction to nonlinear finite element
%              analysis. Oxford University Press, 2004, doi:
%              10.1093/acprof:oso/9780198525295.001.0001
%              [2] S. A. Ragab and H. E. Fayed, Introduction to finite element
%              analysis for engineers. Taylor & Francis Group, 2017
%
%   Part of the PPR repository.
%%
if nargin < 4
    mu = -1;
    if nargin < 3
        lambda = 1;
        if nargin < 2
            eps = 1;
            if nargin < 1
                numElements = 8;
            end
        end
    end
end

%% Prepare a rectangular domain mesh parameters
nnpe=4;             % number of nodes per element (bilinear rectangular element)
a=1;                % length of side in x-direction
b=1;                % length of side in y-direction
nxe=numElements;    % number of elements in x-direction
nye=numElements;    % number of elements in y-direction
nx=nxe+1;           % number of nodes in x-direction
ny=nye+1;           % number of nodes in y-direction
nel=nxe*nye;        % number of elements
nng=nx*ny;          % number of global nodes
nvpn=1;             % number of variables per node
nvpe=nnpe*nvpn;     % number of variables per element
nvg=nng*nvpn;       % number of variables in global mesh

%% Specify material properties and forcing
% Material properties gamma, P(2x2), q(x,y)
gamma = 1;          % mass property that has been nondimensionalized
P = eps*eye(2);     %  ε, diffusion coefficient
q = -lambda;        % -λ, to do with heat generation due to electrical resistance
q3 = -mu;           % -μ, cubic source/reaction term

% Construct forcing function f(x,y) (element-wise constant)
fo = 0; % Currently no force; leaving in case I want to modify it
fgen=zeros(nel,1);
igen=nel;
for ie=1:igen
    fgen(ie)=fo;
end

%% Generate a rectangular mesh on a rectangular domain
% Get connectivity matrix and node locations
[mconn,xyg] = generateRectangularMesh(a,b,nxe,nye,false);

%% Generate element matrices Me, Ke, Fe, and Re and Assemble global system Mg, Kg, Fg, Rg
[Mg,K1g,K2g,K3g,Fg,Rg] = assembleGlobalSystem(gamma,P,q,q3,fgen,mconn,xyg);

%% Prescribe and apply boundary conditions
% Compared to static problems, dynamic problems are a little trickier when
% it comes to imposing boundary conditions. This is because Dirichlet
% boundary conditions reduce the number of unknowns, since the solution
% becomes fixed at those points. They can still be handled, but they are
% just a bit more complicated. So here we will impose Neumann boundary
% conditions on all of the sides to avoid this. The Neumann boundary
% conditions simply prescribe the secondary variable on the RHS, so the
% mass and stiffness are not complicated.
[Mg,K1g,K2g,K3g,Fg,Rg] = applyBoundaryConditions(Mg,K1g,K2g,K3g,Fg,Rg,mconn,xyg,nx,ny,a,b,nel);

%% Put into standard control system form
% The current ODE we have for the global system is
%       Mg u̇ + K1g u + K3g (u ⊗ u ⊗ u) = Fg + Rg.
% We wish to put this in the descriptor state-space form for polynomial
% control systems:
%       Eẋ = A x + F₂ (x ⊗ x) + F₃ (x ⊗ x ⊗ x) + B u
%           y = C x
% Multiplying my M⁻¹ achieves this. The new system will be
%       ̇Mg u = -Kg u + Fg + Rg
%        := A u  + F
% For computational purposes, we will do this inversion by computing the
% Cholesky factor of M; if we were to do M\K, M\F, and M\R these would all
% solve using the Cholesky factor, so we are recycling it. (We could
% simulate faster using sparse matrices and specifying mass matrix as an
% ode option; the inversion destroys the sparsity.)

A = -K1g;
% F2 = -Mchol.' \ (Mchol \ K2g);
F2 = K2g;
% F3 = -K3g;
F3 = K3g; % made negative at element level, saves a lot of time

% Now, changing notation a bit, we want to put the system in the form
%       ̇x = f(x) + g(x) u
% so the old u becomes x and we introduce this new u for controls. f(x) is
% easy, it is just f(x) = A x. g(x) is going to be a B matrix; it is going
% to be computed by taking the F matrix and breaking up into four
% components, one for each boundary. We can use the boundary() function to
% get the node numbers for the boundaries. Note, currently Fg is zero; the
% control input u is really going through the qs variables in Rg, so if you
% add "forcing" fgen later this would have to change.
[AB,BC,CD,DA] = boundary(nx,ny);

B = sparse(nvg,4); % m=4 control inputs, one per boundary
B(AB,1) = Rg(AB);
B(BC,2) = Rg(BC);
B(CD,3) = Rg(CD);
B(DA,4) = Rg(DA);

%% Format outputs
E = sparse(Mg);
f = {A, F2, F3};
g = {B};
h = {1/nvg * ones(1,nvg)}; % average over domain


end

%% Helper functions

function [Mg,K1g,K2g,K3g,Fg,Rg] = applyBoundaryConditions(Mg,K1g,K2g,K3g,Fg,Rg,mconn,xyg,nx,ny,a,b,nel)
% Prescribe and subsequently apply boundary conditions
nnpe=4;             % number of nodes per element (bilinear rectangular element)
nvg = length(Fg);

%% Get boundary node indices
[AB,BC,CD,DA] = boundary(nx,ny);

%% Prescribe boundary conditions using us and qs variables
% These will end up forming the Rg vector, which will be broken into B*u
% for 4 boundary control inputs.
us=zeros(nvg,1); jD = []; % Dirichlet BCs specify u=us (should remain untouched)
qs=zeros(nvg,1); jN = []; % Neumann BCs specify secondary variable

% Prescribe boundary condition on side AB (Neumann)
for i=1:length(AB)
    jN(end+1) = AB(i);
    qs(AB(i)) = qs(AB(i)) + 1; % # is just a placeholder, will be controlled by u
end

% Prescribe boundary condition on side BC (Neumann)
for i=1:length(BC)
    jN(end+1) = BC(i);
    qs(BC(i)) = qs(BC(i)) + 1; % # is just a placeholder, will be controlled by u
end

% Prescribe boundary condition on side CD (Neumann)
for i=1:length(CD)
    jN(end+1) = CD(i);
    qs(CD(i)) = qs(CD(i)) + 1; % # is just a placeholder, will be controlled by u
end

% Prescribe boundary condition on side DA (Neumann)
for i=1:length(DA)
    jN(end+1) = DA(i);
    qs(DA(i)) = qs(DA(i)) + 1; % # is just a placeholder, will be controlled by u
end

%% Impose BCs by modifying Kg, Fg, and Rg
% Apply Dirichlet BCs by replacing those equations with u=us
% ** In this script, there should be none
if(~isempty(jD)) % number of nodes with Essential (Dirichlet) boundary conditions
    for j=1:length(jD) % iterate over rows corresponding to prescribed nodes
        i=jD(j);       % get node index
        K1g(i,:) = 0;   % zero out the row
        K1g(i,i) = 1;   % put a 1 in the coefficient matrix
        Fg(i) = us(i); % put the prescribed value in the rhs vector Fg
    end
end

% Apply Neumann BCs by replacing the secondary variable Re
if(~isempty(jN))
    for ie=1:nel % iterate over elements rather than nodes, since it is the edges that matter
        nodes = mconn(ie,1:nnpe); % extract global node numbers of the nodes of element ie
        % Check if this element has a side (or more) on the global
        % boundary where Neumann boundary condition is given
        [NBC,ib] = checkNBC(nodes,jN);
        if NBC
            xe = xyg(nodes,1); ye = xyg(nodes,2); % extract element node locations
            vne = qs(nodes); % extract vn at the element nodes
            Re = BCN2(vne,ib,xe,ye); % get modified secondary variable Re
            Rg(nodes) = Rg(nodes) + Re;    % assemble global weighted secondary variables vector
        end
    end
end
end

function [Mg,K1g,K2g,K3g,Fg,Rg] = assembleGlobalSystem(gamma,P,q,q3,fgen,mconn,xyg)
% Generate element matrices and assemble the global system
% Note: At this stage no BCs will be applied yet
vec = @(X) X(:);
[nel,nnpe] = size(mconn); % number of elements & number of nodes per element
nng=length(xyg);          % number of global nodes
nvpn=1;             % number of variables per node
nvpe=nnpe*nvpn;     % number of variables per element
nvg=nng*nvpn;       % number of variables in global mesh

%% Assemble cubic global matrix
% Compute element matrices outside of loop = assume same for every element
nodes = mconn(1,1:nnpe);
xe = xyg(nodes,1); ye = xyg(nodes,2);
[Me,K1e,K3e,Fe] = mekefe(gamma,P,q,q3,fgen(1),xe,ye);
Re=zeros(nvpe,1);

% Updated this with the help of chatgpt to assemble more efficiently
% Preallocate cell arrays to store triplet data
Icell  = cell(nel, 1); Jcell  = cell(nel, 1);
VcellM = cell(nel, 1); VcellK = cell(nel, 1);

% Preallocate global force and reaction vectors
Fg = zeros(nvg,1); Rg = zeros(nvg,1);

for ie = 1:nel
    % Extract global node numbers of the nodes of element ie
    % and the global x and y-coordinates of the element nodes
    nodes = mconn(ie,1:nnpe);
    % xe = xyg(nodes,1); ye = xyg(nodes,2);
    
    %% Compute element mass Me, stiffness Ke, forcing Fe, and Re
    % If material property (P and q) and/or force f are element dependent,
    % we need to specify these properties here.
    % [Me,K1e,K3e,Fe] = mekefe(gamma,P,q,q3,fgen(1),xe,ye);
    % Re=zeros(nvpe,1);
    
    % Build triplets from local matrix
    [row_idx, col_idx] = ndgrid(nodes, nodes);
    
    Icell{ie} = row_idx(:); Jcell{ie} = col_idx(:);
    VcellM{ie} = Me(:); VcellK{ie} = K1e(:);
    
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
% Form all zero quadratic global matrix directly
% Since it is empty, sparseIJV is best (it stores basically nothing)
K2g = sparseIJV(nvg, nvg^2);

%% Assemble cubic global matrix
% Updated this with the help of chatgpt to assemble more efficiently
% Preallocate cell arrays to store triplet data
Icell = cell(nel, 1);
Jcell = cell(nel, 1);
Vcell = cell(nel, 1);

for ie = 1:nel
    % [~,~,K3e,~] = mekefe(gamma,P,q,q3,fgen(ie),xe,ye);
    % commented out = assume K3e is constant (it was already computed)
    
    % Local node indices
    nodes = mconn(ie, 1:nnpe);
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
% K3g = sparseCSR(I, J, V, nvg, nvg^3);
K3g = sparseIJV(I, J, full(V), nvg, nvg^3);

end

function [mconn,xyg] = generateRectangularMesh(a,b,nxe,nye,verbose)
% Generate a uniform rectangular mesh for a rectangular domain
%
% Inputs: a - length of side a of the domain
%         b - length of side b of the domain
%       nxe - number of elements in x-direction
%       nye - number of elements in y-direction
%   verbose - optional, plot mesh
%
% Outputs: mconn - connectivity matrix
%            xyg - x-y locations for all of the nodes
%

if nargin < 5; verbose = true; end

nx=nxe+1; ny=nye+1;                     % number of nodes in each direction
nel=nxe*nye;                            % total number of elements
nng=nx*ny;                              % total number of nodes

[x,y] = meshgrid(linspace(0,a,nx),linspace(0,b,ny));

% Construct connectivity matrix for the mesh (maps global nodes to element nodes)
mconn = zeros(nel,4); % Initialize connectivity matrix
for j=1:nye
    jg=j;
    for i=1:nxe
        ig=i;
        kg=(jg-1)*nx+ig;
        ie=(j-1)*nxe+i;
        mconn(ie,1)=kg;
        mconn(ie,2)=kg+1;
        mconn(ie,3)=kg+1+nx;
        mconn(ie,4)=kg+nx;
    end
end

% Construct x-y grid
xyg = zeros(nng,2);
for j=1:ny
    for i=1:nx
        ij=(j-1)*nx+i;
        xyg(ij,:)=[x(j,i), y(j,i)];
    end
end

if verbose % optionally plot the mesh for visualization
    figure(1); hold on; %axis equal;
    
    % Plot grid lines
    plot(x,y); plot(x',y');
    
    % Label node numbers
    for i=1:nng
        text(xyg(i,1),xyg(i,2),num2str(i));
    end
    % Label element numbers
    for i=1:nel
        xp=(xyg(mconn(i,1),1)+xyg(mconn(i,2),1)+xyg(mconn(i,3),1)+xyg(mconn(i,4),1))/4;
        yp=(xyg(mconn(i,1),2)+xyg(mconn(i,2),2)+xyg(mconn(i,3),2)+xyg(mconn(i,4),2))/4;
        text(xp,yp,num2str(i));
    end
end
end

function [AB,BC,CD,DA] = boundary(nx,ny)
% Get the global node indices on the boundary of rectangular domain
% The boundary is composed of line segments AB, BC, CD, DA
%
% Inputs:  nx,ny       - number of nodes in x and y directions
%
% Outputs: AB,BC,CD,DA - vectors mapping to the node numbers on the 4 sides
%

AB = 1:nx;
BC = (1:ny)*nx;
CD = (1:nx) + nx*(ny-1);
DA = (1:ny)*nx - nx + 1;
end

function [Me,K1e,K3e,Fe] = mekefe(gamma,P,q,q3,f,xe,ye)
% Return stiffness matrix and force vector (bilinear rectangular element)
% This function assumes elementwise constant coefficients P, q, and f.
%
% Inputs: gamma - mass property
%             P - material property
%             q - heat convection
%             f - distributed heat source
%         xe,ye - x and y coordinates of element nodes; just for computing
%                 element widths a and b
%
% Outputs:  Me  - element mass matrix
%           K1e - element linear stiffness matrix
%           K3e - element cubic stiffness matrix
%           Fe  - element force vector
%

a=xe(2)-xe(1); b=ye(4)-ye(1);

L11 = b/(6*a)*[ 2 -2 -1  1;
    -2  2  1 -1;
    -1  1  2 -2;
    1 -1 -2  2];
L22 = a/(6*b)*[ 2  1 -1 -2;
    1  2 -2 -1;
    -1 -2  2  1;
    -2 -1  1  2];
L00 = a*b/36*[ 4  2  1  2;
    2  4  2  1;
    1  2  4  2;
    2  1  2  4];

% Cubic stiffness matrix for one element
L3 = a*b/1200 * sparse([48, 36, 9, 36, 0, 24, 12, 18, 0, 0, 4, 12, 0, 0, 0, 24, 0, 0, 0, 0, 0, 12, 9, 6, 0, 0, 6, 8, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12;
    12, 24, 6, 9, 0, 36, 18, 12, 0, 0, 6, 8, 0, 0, 0, 6, 0, 0, 0, 0, 0, 48, 36, 9, 0, 0, 24, 12, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 9, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3;
    3, 6, 4, 6, 0, 9, 12, 8, 0, 0, 9, 12, 0, 0, 0, 9, 0, 0, 0, 0, 0, 12, 24, 6, 0, 0, 36, 18, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 48, 36, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12;
    12, 9, 6, 24, 0, 6, 8, 12, 0, 0, 6, 18, 0, 0, 0, 36, 0, 0, 0, 0, 0, 3, 6, 4, 0, 0, 9, 12, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 24, 0, 0, 0, 36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 48]);


Me = gamma*L00;
K1e = P(1,1)*L11+P(2,2)*L22+q*L00;
Fe = a*b/4 * f*[1; 1; 1; 1];
K3e = q3 * L3;
end

function [NBC,ib] = checkNBC(nodes,jN)
% Check if any of the edges of the element are on the global boundary
% where Natural BCs are prescribed
%
% Inputs: nodes - global node numbers of the current element
%            jN - collection of global node numbers on the boundary
%
% Output: NBC - boolean, true if on a boundary where NBC prescribed
%          ib - vector describing which boundary has NBC prescribed
%
ind1=0; ind2=0; ind3=0; ind4=0;
for k=1:length(jN)
    if(nodes(1)==jN(k)); ind1=1; end
    if(nodes(2)==jN(k)); ind2=1; end
    if(nodes(3)==jN(k)); ind3=1; end
    if(nodes(4)==jN(k)); ind4=1; end
end
ib = [ind1*ind2; ind2*ind3; ind3*ind4; ind4*ind1];

NBC = (sum(ib) > 0); % If any ib not zero, NBC is true
end

function Re = BCN2(vne,ib,xe,ye)
% Compute of Re for a natural boundary condition.
% Assumes linear variation of secondary variable vne on element edges
h12=xe(2)-xe(1);
a=vne(1);
b=vne(2);
q12=h12*[2*a+b; a+2*b; 0; 0];
%%%%%
h23=ye(3)-ye(2);
a=vne(2);
b=vne(3);
q23=h23*[0; 2*a+b; a+2*b; 0];
%%%%%%
h34=xe(3)-xe(4);
a=vne(3);
b=vne(4);
q34=h34*[0; 0; 2*a+b; a+2*b];
%%%%%%%%
h41=ye(4)-ye(1);
a=vne(1);
b=vne(4);
q41=h41*[2*a+b;  0; 0; a+2*b];
%%%%%%%%
Re=(ib(1)*q12+ib(2)*q23+ib(3)*q34+ib(4)*q41)/6;
end


