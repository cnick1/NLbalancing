function x0 = runExample6_getStaticDeflectionIC_sparse(numEls, U0, x0init)
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
    numEls = 16
    U0 = 2e4
    x0init = zeros(6*numEls,1)
end
set(groot,'defaultLineLineWidth',2,'defaultTextInterpreter','latex')
vec = @(x) x(:);
clear n nn F2i F2j F2v F3i F3j F3v I21 I22 I31 I32 I33


n = 6*numEls;
fprintf('Getting Example 6 Initial Condition, n=%d...\n',n)

%% Get dynamics and define control problem
[f,g,h,E] = getSystem6_sparse(numEls,true);

m = size(g{1},2); p = size(h{1},1);
F = @(x) kronPolyEval(f, x);
G = @(x) kronPolyEval(g, x, scenario='G(x)');

% Obtain steady-state Newton iteration for equilibrium point
fsymmetric = f;
for i=2:length(fsymmetric)
    fsymmetric{i} = kronMonomialSymmetrize(fsymmetric{i},n,i);
end

u = [zeros(m-2,1); U0; 0];
if length(x0init) ~= n % interpolate lower-order solution as initial guess
    x0init = vec([zeros(3,2); reshape(x0init,[],2)]); % add fixed node dofs
    x0init = vec(interp1(linspace(0,1,length(x0init)/6), reshape(x0init,[],6), linspace(0,1,numEls+1), 'linear'));
    x0init = reshape(x0init,[],2); % remove fixed node dofs
    x0init = vec(x0init(4:end,:));
end
x0 = newtonIteration(-g{1}*u, @(x) kronPolyEval(f, x), @(x) sparsejcbn(fsymmetric, x), maxIter=10, z0=x0init);
plot(x0); drawnow

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

% k=3 term, need to iterate over n rows to apply kron-vec identity
% xkm1 = kron(x, x);
% for j = 1:n
%         if isempty(find(F{3}(j,:),1)); continue; end
%
%         % J(j,:) = J(j,:) + 3 * xkm1.' * reshape(F{3}(j,:),n^2,[]);
%         J(j,:) = J(j,:) + 3 * (reshape(F{3}(j,:),n^2,[]).' * xkm1).';
% end

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
