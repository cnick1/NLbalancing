function [ft, gt, ht] = transformDynamicsDontInvert(f, g, h, T)
%transformDynamicsDontInvert Transform the model given by f, g, h by transformation T.
%   This function returns the expansions for the transformed dynamics in
%   the form
%             ∂Φ(z)/∂z \dot{z} = ft(z) + gt(z) u,    y = ht(z).
%          NOTE, THE JACOBIAN IS THEN NEEDED ON THE LEFT-HAND-SIDE.
%
%   Usage: [ft, gt, ht] = transformDynamicsDontInvert(f, g, h, T)
%
%   Inputs:
%       f,g,h - cell arrays containing the polynomial coefficients for
%               the drift, input, and output in the original coordinates.
%       T     - cell array containing the polynomial transformation
%               coefficients.
%
%   Output:
%       ft,gt,ht - cell arrays containing the polynomial coefficients
%                  for the transformed drift, input, and output.
%
%   Background: Given a transformation x = Φ(z), we seek to represent the
%    dynamics for the control-affine system
%        \dot{x} = f(x) + g(x) u
%              y = h(x)
%    in the new coordinates as
%        \dot{z} = ftilde(z) + gtilde(z) u
%              y = htilde(z)
%    In general, it is not possible to do this explicitly. Applying the
%    transformation yields
%        ∂Φ(z)/∂z \dot{z} = f(Φ(z)) + g(Φ(z)) u
%                       y = h(Φ(z))
%    and we do not in general have an explicit way to write [∂Φ(z)/∂z]^{-1}.
%    In this function, we will simply return ft(z) = f(Φ(z)), gt(z) =
%    g(Φ(z)), and ht(z) = h(Φ(z)). This can be done exactly for polynomial
%    transformations, though the degree of the dynamics increases, e.g. a
%    cubic function transformed by a quadratic transformation becomes
%    degree 6.
%
%    Note that to simulate the transformed dynamics will require then
%    "inverting the Jacobian" [∂Φ(z)/∂z]^{-1}. The Jacobian matrix can be
%    evaluated using the utility function jcbn(T,z).
%
%   Authors: Nick Corbin, UCSD
%
vec = @(X) X(:);

ld = length(T);
[n,m] = size(g{1}); [p,n] = size(h{1});

lf = length(f); lg = length(g); lh = length(h);
ft = cell(size(f)); gt = cell(size(g)); ht = cell(size(h));

%% Transform drift f(x)
for i = 1:lf*ld
    ft{i} = zeros(n,n^i);
    % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
    % Note: index backwards so that if you only have e.g. T1 you can move on
    for j=flip(1:i)
        if j > lf; continue; end % skip if requesting fj that doesn't exist
        if (i-j)-ld>=0; break; end % skip if requesting T that doesn't exist
        ft{i} = ft{i} + calTTv(T, j, i, f{j}.').';
    end
end

%% Transform input map g(x)
% Instead of dealing with g(x) matrix directly, deal with g_i(x) vectors so
% I can use the same code as for f(x). So loop over m input channels, use gttemp,
% and convert to gt after
gttemp = cell(length(g),m);
for k = 1:m
    for i = 1:(lg-1) % This internal code block is just above for f(x)
        gttemp{i+1,k} = zeros(n,n^i);
        for j=flip(1:i)  % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
            if (i-j)-ld>=0; break; end
            gttemp{i+1,k} = gttemp{i+1,k} + calTTv(T, j, i, g{j+1}(:,k:m:end).').'; % k:m:end needed for converting from g(x)u to sum g_i(x) u_i
        end
    end
end

% Now convert from g_i(x) vectors back to g(x) matrix
gt{1} = g{1}; % B is not state dependent
for i=1:(lg-1)
    gt{i+1} = zeros(n,m*n^i);
    for kk=1:m
        gt{i+1}(:,kk:m:end) = gttemp{i+1,kk};
    end
end

%% Transform output map h(x)
for i = 1:lh
    ht{i} = zeros(p,n^i);
    for j=flip(1:i)  % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
        if (i-j)-ld>=0; break; end
        ht{i} = ht{i} + calTTv(T, j, i, h{j}.').';
    end
end


end
