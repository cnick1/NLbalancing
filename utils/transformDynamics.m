function [ft, gt, ht] = transformDynamics(f, g, h, T)
%transformDynamics Transform the model given by f, g, h by transformation T.
%
%   Usage: [ft, gt, ht] = transformDynamics(f, g, h, T, Tinv)
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
%   Background: Given a transformation T, compute the transformed dynamics.
%   We apply a smooth transformation x = Φ(z) to the control-affine system
%       \dot{x} = f(x) + g(x) u 
%             y = h(x)
%   This yields 
%       ∂Φ(z)/∂z \dot{z} = f(Φ(z)) + g(Φ(z)) u 
%                      y = h(Φ(z))
%    We are seeking the dynamics in the z coordinates, namely 
%       \dot{z} = ftilde(z) + gtilde(z) u 
%    so we expand ftilde(z) and gtilde(z). Inserting these expansions, we
%    then proceed to collect terms of the same degree, leading to the
%    formulas for the coefficients for ftilde(z) & gtilde(z). 
%
%    Half of the challenge is evaluating the transformations f(Φ(z)),
%    g(Φ(z)), and h(Φ(z)). But then we need to "invert the Jacobian", but
%    rather than doing that directly we use the matching of coefficients to
%    directly compute the transformed coefficients. (This is akin to the 
%    projection step where we multiply V' A V). So in practice, the term
%    from inverting the linear transformation coefficient T1 and
%    multiplying by the "transformed" f(Φ(z)), g(Φ(z)) is just the first
%    term we need to consider, but in general there are more. 
%
%   Authors: Nick Corbin, UCSD
%
vec = @(X) X(:);

ld = length(T);
[n,m] = size(g{1}); [p,n] = size(h{1});

lf = length(f); lg = length(g); lh = length(h);
ft = cell(size(f)); gt = cell(size(g)); ht = cell(size(h));

%% Transform drift f(x)
for i = 1:lf
    ft{i} = zeros(n,n^i);
    % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
    % Note: index backwards so that if you only have e.g. T1 you can move on
    for j=flip(1:i) 
        if (i-j)-ld>=0; break; end
        ft{i} = ft{i} + calTTv(T, j, i, f{j}.').';
    end
end

%% Transform input map g(x)
% Instead of dealing with g(x) matrix directly, deal with g_i(x) vectors so 
% I can use the same code as for f(x). So loop over m input channels, use gttemp,
% and convert to gt after
gttemp = cell(length(g),m);
for k = 1:m
    for i = 1:lg-1 % This internal code block is just above for f(x)
        gttemp{i+1,k} = zeros(n,n^i);
        for j=flip(1:i)  % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
            if (i-j)-ld>=0; break; end
            gttemp{i+1,k} = gttemp{i+1,k} + calTTv(T, j, i, g{j+1}(:,k:m:end).').'; % k:m:end needed for converting from g(x)u to sum g_i(x) u_i
        end
    end
end

% Now convert from g_i(x) vectors back to g(x) matrix
gt{1} = g{1}; % B is not state dependent
for i=1:lg-1
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

%% Need to multiply state equation by "inverse of Jacobian", i.e. projection 
% In practice we do not actually multiply by the inverse of the Jacobian,
% we do a coefficient matching that leads to the following sets of terms
% that need to be combined.

% Should also be able to rewrite to use analytical Tinv

% Drift f(x) 
for k = 1:lf
    for i = flip(2:k)                                       % Theoretical sum limits for F_p's; index backwards in i so that p indexes forwards
        p = k + 1 - i;
        if i > ld; continue; end                            % Only run if Ti exists        
        if p > lf; break; end                               % Only run while we have F_p's left
        %%% Naive 
        for ii=1:n 
            ft{k}(ii,:) = ft{k}(ii,:) - LyapProduct(ft{p}.', T{i}(ii,:).', i).';
        end
    end
    ft{k} = T{1}\ft{k};
end

% Input map g(x)
for k = 0:lg-1                                              % 0:lg-1 to deal with zero indexing of Gs
    for i = flip(2:k+1)                                     % Theoretical sum limits for G_p's; index backwards in i so that p indexes forwards
        for jj=1:m
            p = k + 1 - i +1;                                   % Additional +1 to deal with zero indexing of Gs
            if i > ld; continue; end                            % Only run if Ti exists        
            if p > lg; break; end                               % Only run while we have G_p's left
            %%% Naive 
            for ii=1:n
                gt{k+1}(ii,jj:m:end) = gt{k+1}(ii,jj:m:end) - LyapProduct(gt{p}(:,jj:m:end).', T{i}(ii,:).', i).';
            end
        end
    end
    gt{k+1} = T{1}\gt{k+1};
end

% output equation does not get this treatment

end
