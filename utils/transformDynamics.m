function [ft, gt, ht] = transformDynamics(f, g, h, T)
%transformDynamics Transform the model given by f, g, h by transformation T.
%
%   Usage: [ft, gt, ht] = transformDynamics(f, g, h, T)
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
%   TODO: Add more details here.
%
%   Authors: Nick Corbin, UCSD
%
vec = @(X) X(:);

ld = length(T);
[n,m] = size(g{1}); p = length(h{1});

lf = length(f); lg = length(g); lh = length(h);
ft = cell(size(f)); gt = cell(size(g)); ht = cell(size(h));

%% Transform drift f(x)
for i = 1:lf
    ft{i} = zeros(n,n^i);
    for j=1:i % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
        ft{i} = ft{i} + calTTv(T, j, i, f{j}.').';
    end
end

%% Transform input map g(x)
% Instead of dealing with g(x) matrix directly, deal with g_i(x) vectors so 
% I can use the same code as for f(x). So loop over m input channels, use gttemp,
% and convert to gt after
for k = 1:m
    for i = 1:lg-1 % This internal code block is just above for f(x)
        gttemp{i+1,k} = zeros(n,n^i);
        for j=1:i % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
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
    for j=1:i % See Lemma 1 in my SCL paper: Ptilde_i = sum_j^i Pj cT_j,i
        ht{i} = ht{i} + calTTv(T, j, i, h{j}.').';
    end
end

end
