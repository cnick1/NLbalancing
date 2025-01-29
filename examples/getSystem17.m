function [f, g, h] = getSystem17(degree, n)
%getSystem17 Generates a polynomial approximation to an n mass spring
%damper system. The returned model is of size 2*n.
%
%   Usage:  [f,g,h] = getSystem17(degree, n)
%
%   Inputs:
%       degree -  degree d taylor approximation
%       n      -  number of masses
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   References:
%
%%

if nargin < 2
    n = 3;
    if nargin < 1
        degree = 3;
    end
end

filename = fullfile('examples', 'systemData', sprintf('system17_d=%i_n=%i.mat', degree, 2 * n));

if isfile(filename) % load saved system
    load(filename);
else % construct it; can be slow since it uses symbolic calcs
    x = sym('x', [1, 2 * n]).';
    q = x(1:n);
    qdot = x(n + 1:2 * n);
    syms(x);
    
    ms = ones(n, 1); % Warning; if you change this you need to uncomment and use the symbolic code
    ks = ones(n, 1); % Warning; if you change this you need to uncomment and use the symbolic code
    ds = ones(n, 1);
    
    M = diag(ms);
    Minv = diag(1 ./ ms);
    D = diag(ds);
    K = spdiags([-ks, 2 * ks, -ks], [-1 0 1], n, n);
    
    %% Symbolic Lagrangian approach; symbolic computations very slow
    %     T = 0.5 * qdot.' * M * qdot;
    %
    %     %     V = -ks(1)*cos(x(1));
    %     V = -ks(1)*cosP(x(1));
    %     for i=1:n-1
    %         %         V = V - ks(i+1) * cos(x(i)-x(i+1));
    %         V = V - ks(i+1) * cosP(x(i)-x(i+1));
    %     end
    %     %     V = V - ks(n)*cos(x(n)) + n+1; % Ground connection at end
    %     V = V - ks(n)*cosP(x(n)); % Ground connection at end
    %
    %     % simplify(0.5 * q.' * K * q) % Should be linear spring term
    %
    %     fsym = [qdot;
    %         Minv * (gradient(T - V,q) - D * qdot)];
    %     gsym = eye(2*n);
    %     hsym = x;
    %
    %     [f] = approxPolynomialDynamics(fsym,gsym,hsym, x, degree);
    
    %% Direct approach; since we know the analytical Taylor expansion, we code it directly speed improvement
    f{1} = [zeros(n), eye(n); -Minv * K, -D];
    f{2} = sparse(2 * n, (2 * n) ^ 2);
    
    %%
    f{3} = sparse(2 * n, (2 * n) ^ 3);
    
    tensorSize = ones(1, 3) * 2 * n;
    % i=1
    idx = tt_sub2ind(tensorSize, [1 1 1]); % x_1^3/3
    f{3}(1 + n, idx) = 1/3;
    idx = tt_sub2ind(tensorSize, [2 2 2]); % -x_2^3/6
    f{3}(1 + n, idx) = -1/6;
    idx = tt_sub2ind(tensorSize, [2 2 1]); % (x1*x2^2)/2
    f{3}(1 + n, idx) = 1/2;
    idx = tt_sub2ind(tensorSize, [2 1 1]); % -(x1^2*x2)/2
    f{3}(1 + n, idx) = -1/2;
    
    for i = 2:n - 1
        idx = tt_sub2ind(tensorSize, [i + 1 i + 1 i + 1]); % -x_{i+1}^3/6
        f{3}(i + n, idx) = -1/6;
        idx = tt_sub2ind(tensorSize, [i + 1 i + 1 i]); % (x_i*x_{i+1}^2)/2
        f{3}(i + n, idx) = 1/2;
        idx = tt_sub2ind(tensorSize, [i + 1 i i]); % -(x_i^2*x_{i+1})/2
        f{3}(i + n, idx) = -1/2;
        idx = tt_sub2ind(tensorSize, [i i i]); % x_i^3/3
        f{3}(i + n, idx) = 1/3;
        idx = tt_sub2ind(tensorSize, [i i i - 1]); % -(x_{i-1}*x_i^2)/2
        f{3}(i + n, idx) = -1/2;
        idx = tt_sub2ind(tensorSize, [i i - 1 i - 1]); % (x_{i-1}^2*x_i)/2
        f{3}(i + n, idx) = 1/2;
        idx = tt_sub2ind(tensorSize, [i - 1 i - 1 i - 1]); % - x_{i-1}^3/6
        f{3}(i + n, idx) = -1/6;
    end
    
    % i=n
    idx = tt_sub2ind(tensorSize, [n n n]); % x_n^3/3
    f{3}(2 * n, idx) = 1/3;
    idx = tt_sub2ind(tensorSize, [n - 1 n - 1 n - 1]); % -x_n-1^3/6
    f{3}(2 * n, idx) = -1/6;
    idx = tt_sub2ind(tensorSize, [n n - 1 n - 1]); % (x1*x2^2)/2
    f{3}(2 * n, idx) = 1/2;
    idx = tt_sub2ind(tensorSize, [n n n - 1]); % -(x1^2*x2)/2
    f{3}(2 * n, idx) = -1/2;
    
    %     g = speye(2 * n); h = speye(2 * n);
    %     g = diag([zeros(n,1);ones(n,1)]); h = diag([ones(n,1);zeros(n,1)]);
    g = spdiags([zeros(n,1);ones(n,1)],0,2*n,2*n); h = spdiags([ones(n,1);zeros(n,1)],0,2*n,2*n);
    
    %% Save system data so it doesn't need to be computed every time
    save(filename, 'f', 'g', 'h')
end
end

% function cosPS = cosP(y)
% cosPS = - y^2/2 + y^4/24;
% end
