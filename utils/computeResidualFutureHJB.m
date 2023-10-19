function RES = computeResidualFutureHJB(f, g, h, eta, w, degree, dataRange, nData)
if nargin < 8
    nData = 301;
    if nargin < 7
        dataRange = 1;
        if nargin < 6
            degree = length(w);
        end
    end
end
w = w(1:degree);

n = length(f{1}); % Get state dimension from A matrix

if n > 1 % Remove extra zeros if the model has been tweaked; these should all be matrices
    f = f(~cellfun(@isscalar, f));
    g = g(~cellfun(@isscalar, g));
    h = h(~cellfun(@isscalar, h));
end

% Construct inf-norm unit ball, i.e. set up x,y,z,... from (-1,1) w/ N pts
N = nData; % Number of values between -1 and 1

% Generate linspace(-1, 1, N) once
xn = linspace(-dataRange, dataRange, N);

% Iterate through all points in the discretized space and compute the RES
RES = zeros(N ^ n, 1);
for i = 1:N ^ n
    % Calculate the indices for each dimension
    indices = mod(floor((i - 1) ./ N .^ (0:(n - 1))), N) + 1;
    x = flip(xn(indices).'); % This is the ith point x in the state-space

    if length(g) > 1
        % Polynomial input
        RES(i) = (0.5 * kronPolyDerivEval(w, x)) * kronPolyEval(f, x) ...
            - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * (g{1} + kronPolyEval(g(2:end), x)) * (g{1} + kronPolyEval(g(2:end), x)).' * kronPolyDerivEval(w, x).' ...
            + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);
    else
        % Linear/constant input B
        RES(i) = (0.5 * kronPolyDerivEval(w, x)) * kronPolyEval(f, x) ...
            - eta / 2 * 0.25 * kronPolyDerivEval(w, x) * g{1} * g{1}.' * kronPolyDerivEval(w, x).' ...
            + 0.5 * kronPolyEval(h, x).' * kronPolyEval(h, x);
    end

end

if n > 1
    RES = reshape(RES, N * ones(1, n));
end
end
