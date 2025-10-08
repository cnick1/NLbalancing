function [f, g, Q, D2, y] = getSystem9Neumann(eps, n)
%getSystem9Neumann Cubic Allen-Cahn w/ Nuemann BCs model using Chebychev discretization [1,2].
%
%   Usage:  [f,g,Q,D,y] = getSystem9Neumann(eps, N, y0)
%
%   Background: The Allen-Cahn equation is
%                   uₜ = eps*uₓₓ+u-u³, uₓ(-1)=0, uₓ(1)=0
%       The spatial domain is discretized using Chebychev points, and the
%       spatial derivative becomes uₓ = D*u with differentiation matrix D.
%       The spatial domain x and differentiation matrix D are given by the
%       cheb() function from [1]. The system has equilibrium solutions
%       u(x) = -1, 0, 1, which all satisfy the Neumann boundary conditions.
%       The family of steady-state solutions given by u(x) = tanh((x-x0)/sqrt(2*eps))
%       also form equilibrium solutions.
%
%       Upon discretization, the model is
%                   uₜ = eps*D²*u+u-u³
%
%           -> A = eps*D² + I
%       The equilibrium at the origin is the unstable u(x)=0, so we seek to
%       stabilize it.
%
%       The model and code is based on p34.m from [1].
%
%   References: [1] L. N. Trefethen, Spectral methods in MATLAB. Society
%               for Industrial and Applied Mathematics, 2000.
%               doi: 10.1137/1.9780898719598.
%               [2] S. Dolgov, D. Kalise, and K. Kunisch, "Tensor
%               decomposition methods for high-dimensional
%               Hamilton-Jacobi-Bellman equations," SIAM Journal on
%               Scientific Computing, vol. 43, no. 3, pp. A1625–A1650, Jan.
%               2021, doi: 10.1137/19m1305136.
%
%   Part of the PPR repository.
%%

if nargin < 3
    y0 = 0.5;
    if nargin < 2
        n = 14;
        if nargin < 1
            eps = 0.01;
        end
    end
end
N=n+1;

% Differentiation matrix
[D,y] = cheb(N); D2 = D^2;
s_bound = -[D(1,1), D(1,N+1); D(N+1,1), D(N+1,N+1)] \ D([1,N+1], 2:N);
bound = [1; N+1];
D2 = D2(2:N, 2:N) + D2(2:N, bound)*s_bound;

f{1} = eps*D2 + eye(n);
f{2} = sparse(1:n,linspace(1,n^2,n),0);
f{3} = sparse(1:n,linspace(1,n^3,n),-1);

% Form B matrix using code from TT-HJB, B was called gxi in the original script
% Boost up approx order by 1 by differentiating an exact integral
omega_left = -0.5;omega_right = 0.2;
B = double(y>=omega_left & y<=omega_right).*(y-omega_left) + double(y>omega_right).*(omega_right-omega_left);
B = D*B;
B = B(2:N);

g = B;
y = y(2:N);

%% Cost state penalty matrix x.'Qx
G0 = D(1:N,1:N)';
wxi = G0\eye(N,1);
wxi = wxi(2:N) + wxi(1)*s_bound(1,:)';
Q = diag(wxi);

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