function w = approxFutureEnergy(f, g, h, eta, d, verbose)
%approxFutureEnergy  Compute a polynomial approximation to the future
% energy function for a polynomial control-affine dynamical system.
%
%   Usage: w = approxFutureEnergy(f,g,h,eta,d,verbose)
%
%   Inputs:
%       f,g,h   - cell arrays containing the polynomial coefficients
%                 for the drift, input, and output.
%                   • f must contain at least linear and quadratic coefficients
%                   • g must contain at least a linear input (B matrix)
%                   • h must contain at least a linear input (C matrix)
%       eta     - η=1-1/γ^2, where γ is the H∞ gain parameter. For open-loop
%                 balancing, use eta=0. For closed-loop (HJB) balancing, use
%                 eta=1. Any other value between -1 and ∞ corresponds to
%                 H∞ balancing.
%       d       - desired degree of the computed energy function. A degree d
%                 energy function uses information from f,g,h up-to degree d-1.
%                 The default choice of d is lf+1, where lf is the degree of
%                 the drift.
%       verbose - optional argument to print runtime information
%
%   Output:
%       w       - cell array containing the polynomial energy function coefficients
%
%   Background: Computes a degree d polynomial approximation to the energy function
%
%          E^+(x) = 1/2 ( w{2}'*(x⊗x) + ... + w{d}'*(...⊗x) )
%
%   for the polynomial control-affine system
%
%    \dot{x} = Ax + F2*(x⊗x) + F3*(x⊗x⊗x) + ...
%              + Bu + G1*(x⊗u) + G2*(x⊗x⊗u) + ...
%          y = Cx + H2*(x⊗x) + H3*(x⊗x⊗x) + ...
%
%   where eta = η=1-1/γ^2, where γ is the H∞ gain parameter. w{2} = vec(W2) = W2(:)
%   solves the H∞ Algebraic Riccati Equation
%
%    A'*W2 + W2*A - eta*W2*B*B'*W2 + C'*C = 0,
%
%   and the remaining w{i} solve linear systems arising from the Future H∞
%   Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%   Details are in Section III.B of reference [1] or III.A of reference [2].
%
%   Requires the following functions from the KroneckerTools repository:
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%
%   Authors: Jeff Borggaard, Virginia Tech
%            Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear
%               balanced truncation: Part 1—computing energy functions,” arXiv,
%               Dec. 2022. doi: 10.48550/ARXIV.2209.07645
%              [2] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞
%               energy functions for polynomial control-affine systems,” 2023.
%
%             See Algorithm 1 in [1].
%
%  Part of the NLbalancing repository.
%%

if (nargin < 6)
    verbose = false;
end

% Create pointer/shortcuts for dynamical system polynomial coefficients
if iscell(f)
    % polynomial drift balancing
    A = f{1};
    N = f{2}; % maybe don't do this here? Well if N is missing the code would break anyways
    lf = length(f);
else
    error("Must pass in at least quadratic dynamics")
end

if iscell(g)
    % QB or polynomial input balancing
    B = g{1};
    lg = length(g) - 1;
else
    % Will reduce to Jeff's original code
    B = g;
    lg = 0;
    g = {B};
end

if iscell(h)
    % polynomial output balancing
    C = h{1};
    lh = length(h);
else
    % Will reduce to Jeff's original code
    C = h;
    lh = 1;
    h = {C};
end

n = size(A, 1); % A should be n-by-n
m = size(B, 2); % B should be n-by-m
p = size(C, 1); % C should be p-by-n

% Create a vec function for readability
vec = @(X) X(:);

%% k=2 case
R = eye(m) / eta;

if (eta > 0)
    [W2] = icare(A, B, (C.' * C), R);

    if (isempty(W2) && verbose)
        warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
    end

elseif (eta < 0)
    [W2] = icare(A, B, (C.' * C), R, 'anti');

    if (isempty(W2) && verbose)
        warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
        warning('approxFutureEnergy: using the hamiltonian')
        [~, W2, ~] = hamiltonian(A, B, C.' * C, R, true);
    end

else % eta==0
    [W2] = lyap(A.', (C.' * C));

end

if (isempty(W2))
    error('approxFutureEnergy: Can''t find a stabilizing solution')
end

%  Check the residual of the Riccati/Lyapunov equation
if (verbose)
    RES = A' * W2 + W2 * A - eta * (W2 * B) * (B' * W2) + C' * C;
    fprintf('The residual of the Riccati equation is %g\n', norm(RES, 'inf'));
    clear RES
end

%  Reshape the resulting quadratic coefficients
w{2} = vec(W2);

%% k=3 case
if (d > 2)
    GaWb = cell(2 * lg + 1, d - 1); % Pre-compute G_a.'*W_b, etc for all the a,b we need
    GaWb{1, 2} = B.' * W2;
    % set up the generalized Lyapunov solver
    [Acell{1:d}] = deal(A.' - eta * W2 * (B * B.'));

    b = -LyapProduct(N.', w{2}, 2);

    if lg > 0 % New for QB/polynomial input
        Im = speye(m);
        GaWb{2, 2} = g{2}.' * W2;
        b = b + 2 * eta * kron(speye(n ^ 3), vec(Im).') * vec(kron(GaWb{2, 2}, GaWb{1, 2}));
    end

    % New for polynomial output h(x)
    % TODO: use symmetry to cut in half
    for p = (3 - lh):lh
        q = 3 - p;
        b = b - vec(h{p}.' * h{q}');
    end

    [w{3}] = KroneckerSumSolver(Acell(1:3), b, 3); % Solve Ax=b for k=3
    [w{3}] = kronMonomialSymmetrize(w{3}, n, 3); % Symmetrize w3

    %% k>3 cases (up to d)
    for k = 4:d
        GaWb{1, k - 1} = B.' * reshape(w{k - 1}, n, n ^ (k - 2));

        b = -LyapProduct(N.', w{k - 1}, k - 1); % Pre-compue all the L(N') terms

        % New for polynomial drift f(x)

        iRange = 2:(k - 2);
        iRange = iRange(max(k - lf, 1):end); % Need to only do lf last i's; if there are only 2 Ns for example, only need k-2! otherwise f(xi) doesnt exist and would require a bunch of empty f(xi)s
        %         TODO: may be better to write directly in terms of xi, but then how to get rid of i=k-1...?
        for i = iRange % would be from 2:k-1 but k-1 was covered in instantiation of b
            xi = k + 1 - i;
            b = b - LyapProduct(f{xi}.', w{i}, i);
        end

        % Now add all the terms from the 'B' sum by looping through the i and j
        for i = 3:(k + 1) / 2 % i+j=k+2
            j = k + 2 - i;
            tmp = GaWb{1, i}.' * GaWb{1, j};
            b = b + 0.25 * eta * i * j * (vec(tmp) + vec(tmp.'));
        end

        if ~mod(k, 2) % k is even
            i = (k + 2) / 2;
            j = i;
            tmp = GaWb{1, i}.' * GaWb{1, j};
            b = b + 0.25 * eta * i * j * vec(tmp);
        end

        % Now add the higher order polynomial terms "G" by iterating through the sums
        [g{lg + 2:2 * lg + 1}] = deal(0); % Need extra space in g because of GaWb indexing

        for o = 1:2 * lg
            for idx = 2:k - 1 % Might be repetitive
                if o + 1 < lg + 2
                    GaWb{o + 1, idx} = g{o + 1}.' * sparse(reshape(w{idx}, n, n ^ (idx - 1)));
                else
                    GaWb{o + 1, idx} = 0;
                end
            end
            for p = max(0, o - lg):min(o, lg)

                for i = 2:k - o
                    q = o - p;
                    j = k - o - i + 2;
                    tmp = kron(speye(n ^ p), vec(Im).') ...
                        * kron(vec(GaWb{q + 1, j}).', kron(GaWb{p + 1, i}, Im)) ...
                        * kron(speye(n ^ (j - 1)), kron(perfectShuffle(n ^ (i - 1), n ^ q * m), Im)) ...
                        * kron(speye(n ^ (k - p)), vec(Im));
                    b = b + 0.25 * eta * i * j * vec(tmp);
                end

            end

        end

        % New for polynomial output h(x)
        % TODO: use symmetry to cut in half
        for p = (k - lh):lh % would be 1:(k-1) but need to truncate only to h{} terms which exist
            q = k - p;
            b = b - vec(h{p}.' * h{q}');
        end

        % Done with RHS! Now solve and symmetrize!
        [w{k}] = KroneckerSumSolver(Acell(1:k), b, k);
        [w{k}] = kronMonomialSymmetrize(w{k}, n, k);
    end

end

% if verbose % Check HJB Residual
%     % Compute "residual":
%     %     nX = 301; nY = nX;
%     %     xPlot = linspace(-1, 1, nX);
%     %     yPlot = linspace(-1, 1, nY);
%     %     [X, Y] = meshgrid(xPlot, yPlot);
%     %     RES = zeros(nY, nX);
%     %     degree = length(w);
%     %     for i = 1:nY
%     %         for j = 1:nX
%     %             x = [X(i, j); Y(i, j)];
%     %
%     %             RES(i, j) = (0.5 * kronPolyDerivEval(w, x, degree)) * (f{1} * x + f{2} * kron(x, x)) ...
%     %                 - eta / 2 * 0.25 * kronPolyDerivEval(w, x, degree) * (g{1} + g{2} * x) * (g{1} + g{2} * x).' * kronPolyDerivEval(w, x, degree).' ...
%     %                 + 0.5 * (C * x).' * (C * x);
%     %         end
%     %     end
%     % RES = sum(sum(abs(RES))) / (nX * nY);
%     RES = computeResidualFutureHJB(f, g, h, eta, w);
%
%     fprintf('The residual of the HJB equation on the unit square is %g\n', norm(RES, 'inf'));
% else
%     RES = [];
% end

end
