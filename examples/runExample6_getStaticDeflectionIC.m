function x0s = runExample6_getStaticDeflectionIC(numEls, U0)
%runExample6_getStaticDeflectionIC
%
%   Usage:  runExample6_getStaticDeflectionIC(degree,lim)
%
%   Inputs:    degree - desired degree of the energy function approximation
%                 lim - the size of the grid in the z coordinates
%
%   Description:
%
%   References:
%
%   Part of the NLbalancing repository.
%%
% close all;
arguments
    numEls = [1 2 4 8 16 32 64 128 180]
    U0 = 2e4
end

set(groot,'defaultLineLineWidth',2,'defaultTextInterpreter','latex')
IC_filename = sprintf('examples/getSystem6_staticDeflectionIC_%i.mat',U0);

if isfile(IC_filename)
    load(IC_filename)
else
    x0s = cell(length(numEls),1);
    x0s{1} = getX0(numEls(1),U0,zeros(6*numEls(1),1));
    for i = 2:length(numEls)
        x0s{i} = getX0(numEls(i), U0, x0s{i-1});
    end
    save(IC_filename, 'x0s');
end

end

function x0 = getX0(numEls,U0,x0init)
vec = @(x) x(:);
clear n nn F2i F2j F2v F3i F3j F3v I21 I22 I31 I32 I33

n = 6*numEls;
fprintf('Getting Example 6 Initial Condition, n=%d...\n',n)

%% Get dynamics and define control problem
[f,g,h,E] = getSystem6_sparse(numEls,true);

m = size(g{1},2); p = size(h{1},1);
F = @(x) kronPolyEval(f, x);
G = @(x) g{1};

% Obtain steady-state Newton iteration for equilibrium point
fsymmetric = f;
for i=2:length(fsymmetric)
    fsymmetric{i} = sparseKronMonomialSymmetrize(fsymmetric{i},n,i);
end

u = [zeros(m-2,1); U0; 0];
if length(x0init) ~= n % interpolate lower-order solution as initial guess
    x0init = vec([zeros(3,2); reshape(x0init,[],2)]); % add fixed node dofs
    x0init = vec(interp1(linspace(0,1,length(x0init)/6), reshape(x0init,[],6), linspace(0,1,numEls+1), 'linear'));
    x0init = reshape(x0init,[],2); % remove fixed node dofs
    x0init = vec(x0init(4:end,:));
end
x0 = newtonIteration(-g{1}*u, @(x) kronPolyEval(f, x), @(x) sparsejcbn(fsymmetric, x), maxIter=10, z0=x0init);

end


function [c] = sparseKronMonomialSymmetrize(c, n, k)
if nargin < 3
    k = log(length(c))/log(n);
end

if ~any(c(:))
    return
end


%% Perform actual symmetrization
for i=1:size(c,1)
    if any(c(i,:))
        I = find(c(i,:));
        subs = tt_ind2sub(ones(1, k) * n, I);

        all_perms = cellfun(@perms, num2cell(unique(subs, 'rows'), 2), 'UniformOutput', false);
        subs = unique(vertcat(all_perms{:}), 'rows');
        idx = tt_sub2ind(ones(1, k) * n, subs);

        classsubs = sort(subs, 2); % Normalize to one permutation, i.e. reference element

        mult = [1 cumprod(ones(1, k - 1) * n)]; % Form shifts
        linclassidx = (classsubs - 1) * mult' + 1; % Form vector that maps to the reference elements

        classsum = sparse(linclassidx, 1, c(i, idx).', n^k, 1);
        classnum = sparse(1:max(linclassidx),1,accumarray(linclassidx, 1),n^k,1);

        [I,J,~] = find(classsum);
        avg = sparse(I,J,classsum(I) ./ classnum(I),n^k,1);

        % Copy averaged entries back to the other locations
        avg(idx) = avg(linclassidx);

        % Fill in each entry with its new symmetric version
        c(i,:) = avg;
    end
end
end


function v = symmetrizeHelper(v,linclassidx,classnum)
%symmetrizeHelper Performs the gather/scatter averaging operation
%   Input/Output: v - vector of dimension nᵏ, symmetrized result
if ~any(v)
    return
end

classsum = accumarray(linclassidx, v);

% Option 2:
if issparse(v)
    [I,J,~] = find(classsum);
    avg = sparse(I,J,classsum(I) ./ classnum(I),length(v),1);
else
    avg = zeros(size(v));
    I = find(classsum);
    avg(I) = classsum(I) ./ classnum(I);
end


% Fill in each entry with its new symmetric version
v = avg(linclassidx);
end

function J = sparsejcbn(F, x)
%jcbn Return the Jacobian J(x) = ∂f(x)/∂x of the function f(x) evaluated at x.
%              f(x) = F₁x + F₂(x⊗x) + ... + Fd(x...⊗x)
%   the Jacobian is given by
%       J(x) = ∂f(x)/∂x = F₁ + 2F₂(I⊗x) + ... + d Fd(I...⊗x)
%%
n = size(x, 1);
% k=1 term
J = full(double(F{1}));
if isempty(find(x,1)) % if x is zero, only the constant term remains, so we are done
    return;
end
% k=2 term, need to iterate over n rows to apply kron-vec identity
for j = 1:n
    if isempty(find(F{2}(j,:),1)); continue; end
    J(j,:) = J(j,:) + 2 * x.' * reshape(F{2}(j,:),n,[]);
end

% Option 2:
persistent nn F3i F3j F3v I1 I2
if isempty(nn) || nn ~= length(x)
    nn = length(x);
    for j = 1:nn
        if isempty(find(F{3}(j,:),1)); continue; end
        [F3i{j}, F3j{j}, F3v{j}] = find(reshape(F{3}(j,:),nn,[]));
        [I1{j}, I2{j}] = ind2sub([nn nn], F3j{j});
    end
    
end
for j = 1:n
    % if isempty(find(F{3}(j,:),1)); continue; end
    if isempty(find(F{3}(j,:),1)); continue; end
    xprod = x(I1{j}) .* x(I2{j});
    J(j,:) = J(j,:) + 3 * accumarray(F3i{j}, F3v{j} .* xprod, [n, 1]).';
end

end
