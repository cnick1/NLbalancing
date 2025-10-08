function [f, g, h, xyg] = getSystem28(numElements)
%getSystem28  Generates a linear 2D cubic finite element heat equation model system.
%            The system is a finite element model for a LINEAR heat equation.
%            The function returns a finite element model with bilinear
%            rectangular elements. The model is subject to boundary control
%            via the secondary variables, so there are 4 inputs. The output
%            is the mean temperature. Goal is to bring the system to zero
%            from some initial condition. 
%
%   Usage:   [f,g,h] = getSystem28()
%
%   Inputs:
%       numElements    - number of elements to discretize each direction with
%                        (default = 32, leading to n=1089 dimension model)
%
%   Outputs:    f,g,h  - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%                  xyg - global node locations
%
%   Description: after finite element discretization, the finite element
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
%   Reference: [1] J. N. Reddy, An introduction to nonlinear finite element
%              analysis. Oxford University Press, 2004,
%              doi: 10.1093/acprof:oso/9780198525295.001.0001
%              [2] S. A. Ragab and H. E. Fayed, Introduction to finite
%              element analysis for engineers. Taylor & Francis Group, 2017
%
%   Part of the NLbalancing repository.
%%

vec = @(X) X(:);

if nargin < 1
    numElements = 32;
end

%% Prepare a rectangular domain mesh parameters
nnpe=4;             % number of nodes per element (bilinear rectangular element)
a=1;              % length of side in x-direction
b=1;              % length of side in y-direction
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
gamma = 41.9404;
P = 5.03e-3 * eye(2); % thermal conductivity matrix, k_xx & k_yy
q = .1;        % some heat convection for a stable steady state

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

%% Generate element matrices Ke, Fe, and Re and Assemble global system Kg, Fg, Rg
[Mg,Kg,Fg,Rg] = assembleGlobalSystem(gamma,P,q,fgen,mconn,xyg);

%% Prescribe and apply boundary conditions
% Compared to static problems, dynamic problems are a little trickier when
% it comes to imposing boundary conditions. This is because Dirichlet
% boundary conditions reduce the number of unknowns, since the solution
% becomes fixed at those points. They can still be handled, but they are
% just a bit more complicated. So here we will impose Neumann boundary
% conditions on all of the sides to avoid this. The Neumann boundary
% conditions simply prescribe the secondary variable on the RHS, so the
% mass and stiffness are not complicated. 
[Mg,Kg,Fg,Rg] = applyBoundaryConditions(Mg,Kg,Fg,Rg,mconn,xyg,nx,ny,a,b,nel);

%% Put into standard control system form 
% The current ODE we have for the global system is Mg u̇ + Kg u = Fg + Rg. 
% We wish to put this in the standard state-space form for control systems:
%       ẋ = A x + F₂ (x ⊗ x) + F₃ (x ⊗ x ⊗ x) + ...
%               + B u + G₁ (x ⊗ u) + G₂ (x ⊗ x ⊗ u) + ...
%           y = C x + H₂ (x ⊗ x) + ...
% Multiplying my M⁻¹ achieves this. The new system will be 
%       ̇u = -Mg⁻¹Kg u + Mg⁻¹Fg + Mg⁻¹Rg
%        := A u  + F
% For computational purposes, we will do this inversion by computing the 
% Cholesky factor of M; if we were to do M\K, M\F, and M\R these would all 
% solve using the Cholesky factor, so we are recycling it. (We could  
% simulate faster using sparse matrices and specifying mass matrix as an 
% ode option; the inversion destroys the sparsity.)

Mchol = chol(Mg).'; % Use Cholesky factor for inverting 
A = -Mchol.' \ (Mchol \ Kg);
% F =  Mchol.' \ (Mchol \ (Fg + Rg));

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

B = zeros(nvg,4); % m=4 control inputs, one per boundary
B(AB,1) = Rg(AB);
B(BC,2) = Rg(BC);
B(CD,3) = Rg(CD);
B(DA,4) = Rg(DA); 

B =  Mchol.' \ (Mchol \ B);


%% Placeholder for nonlinear terms
F2 = sparse(nvg,nvg^2);
F3 = sparse(nvg,nvg^3);

%% Format outputs
f = {A, F2, F3};
g = {B};
h = {1/nvg * ones(1,nvg)}; % average over domain


end

%% Helper functions

function [Mg,Kg,Fg,Rg] = applyBoundaryConditions(Mg,Kg,Fg,Rg,mconn,xyg,nx,ny,a,b,nel)
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
        Kg(i,:) = 0;   % zero out the row
        Kg(i,i) = 1;   % put a 1 in the coefficient matrix
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

function [Mg,Kg,Fg,Rg] = assembleGlobalSystem(gamma,P,q,fgen,mconn,xyg)
% Generate element matrices and assemble the global system
% Note: At this stage no BCs will be applied yet
[nel,nnpe] = size(mconn); % number of elements & number of nodes per element
nng=length(xyg);          % number of global nodes
nvpn=1;             % number of variables per node
nvpe=nnpe*nvpn;     % number of variables per element
nvg=nng*nvpn;       % number of variables in global mesh

% Initilize matrices
Mg=zeros(nvg,nvg); % global mass matrix
Kg=zeros(nvg,nvg); % global stiffness matrix
Fg=zeros(nvg,1);   % global force vector
Rg=zeros(nvg,1);   % global weighted secondary variable vector

for ie=1:nel
    % Extract global node numbers of the nodes of element ie
    % and the global x and y-coordinates of the element nodes
    nodes = mconn(ie,1:nnpe);
    xe = xyg(nodes,1); ye = xyg(nodes,2);

    %% Compute element mass Me, stiffness Ke, forcing Fe, and Re
    % If material property (P and q) and/or force f are element dependent,
    % we need to specify these properties here.
    [Me,Ke,Fe] = mekefe(gamma,P,q,fgen(ie),xe,ye);
    Re=zeros(nvpe,1);

    % Assemble element matrices Ke, fe, Re into global matrix Kg, Fg, Rg
    Fg(nodes) = Fg(nodes) + Fe;
    Rg(nodes) = Rg(nodes) + Re;
    Kg(nodes,nodes) = Kg(nodes,nodes) + Ke;
    Mg(nodes,nodes) = Mg(nodes,nodes) + Me;
end
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

[y,x] = meshgrid(linspace(0,a,nx),linspace(0,b,ny));

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
        xyg(ij,:)=[x(i,j), y(i,j)];
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

function [Me,Ke,Fe] = mekefe(gamma,P,q,f,xe,ye)
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
% Outputs:   Me - element mass matrix
%            Ke - element stiffness matrix         
%            Fe - element force vector
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

Me = gamma*L00;
Ke = P(1,1)*L11+P(2,2)*L22+q*L00;
Fe = a*b/4 * f*[1; 1; 1; 1];
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


