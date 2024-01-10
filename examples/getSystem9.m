function [f, g, h, D, y] = getSystem9(eps, N, y0)
%getSystem9 Returns a cubic Allen-Cahn model using Chebychev spatial discretization [1].
%
%   Usage:  [f,g,h,D,y] = getSystem9(eps, N, y0)
%
%   Background: The Allen-Cahn equation is
%                   u_t = eps*u_xx+u-u^3, u(-1)=-1, u(1)=1
%       The spatial domain is discretized using Chebychev points, and the
%       spatial derivative becomes u_x = D*u with differentiation matrix D.
%       The spatial domain x and differentiation matrix D are given by the
%       cheb() function from [1]. The system has equilibrium solutions
%       u(x) = -1, 0, 1; however, these don't satisfy the boundary
%       conditions. Another equilibrium solution exists, which gives the
%       steady-state of the system: u(x) = tanh((x-x0)/sqrt(2*eps)).
%
%       Upon discretization, the model is
%                   u_t = eps*D^2*u+u-u^3
%
%           -> A = eps*D^2 + I
%       But the equilibrium at the origin is the unstable u(x)=0. We shift
%       the equilibrium u_ref = tanh((x-x0)/sqrt(2*eps)) to the origin,
%       giving a new A matrix and introducing a quadratic drift component
%       F2 as well. u^3 can be written as F3*kron(u,kron(u,u)).
%
%       The model and code is based on p34.m from [1].
%
%   References: [1] L. N. Trefethen, Spectral methods in MATLAB. Society
%               for Industrial and Applied Mathematics, 2000.
%               doi: 10.1137/1.9780898719598.
%
%%
if nargin < 3
    y0 = [];
    if nargin < 2
        N = 20;
        if nargin < 1
            eps = 0.01;
        end
    end
end

% Differentiation matrix
[D,y] = cheb(N); D2 = D^2;
D2([1 N+1],:) = zeros(2,N+1);

% Construct shifted dynamics with u(x)=tanh((x-x0)/sqrt(2*eps)) equilibrium @ origin
n = N+1;

% Shift reference equilibrium to the origin
% vref = tanh((y-y0)/sqrt(2*eps));
if isempty(y0)
    vref= 0;
else
    load(fullfile('examples', 'systemData',sprintf('system9_equilibria_N=%i.mat',N)))
    vref = d(y0); vref = vref{1}; %plot(y,vref)
end

f{1} = eps*D2 + eye(n) - 3*diag(vref.^2);
f{2} = sparse(1:n,linspace(1,n^2,n),-3*vref);
f{3} = sparse(1:n,linspace(1,n^3,n),-1);
g = eye(n); 
h = eye(n);

% Select controls |--*--*--*--|

end

function [D,x] = cheb(N)
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries
end