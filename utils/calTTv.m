function [Tv] = calTTv(T, m, k, v)
%calTTv Calculates the term \cal{T}_{m,k}.'v
%
%  Note the additional T in calTTv ensures the transpose of T is taken.
%
%         Tv = \|i\|=k kron(T_i1.',kron(T_i2.',kron(T_i3.',...,T_im.'))*v
%
%  Usage:  Tv = calTTv(T,m,k,v)
%
%  Variables:   T  a cell array of matrices.  The pth entry has size n times n^p
%               v  a vector of dimension (n^m)
%
%  This is multiplication performed recursively using Kronecker product rules.
%
%  This function assumes the functions mkIndices and kroneckerRight have
%  been imported from the KroneckerTools repository (and in the Matlab path).
%
%  Author: Jeff Borggaard, VT
%          Nick Corbin, UCSD
%
%  License: MIT
%
%  Part of the NLbalancing repository.
%%
method = 3;

switch method
    case 1 %% Old way
        % Get a list of indices
        indexSet = mkIndices(m, k);
        nTerms = size(indexSet, 1);
        
        n = size(T{1}, 2);
        
        Tv = zeros(1, n ^ k);
        
        for i = 1:nTerms
             Tv = Tv + kroneckerRight(v.', T(indexSet(i, :)));
        end
        
        Tv = Tv.';
        
    case 2 %% basically old way but using new function and perms (just a test)
        % Get a list of indices
        [indexSet, mult] = findCombinations(m, k);
        nTerms = size(indexSet, 1);
        %
        n = size(T{1}, 2);
        
        Tv = zeros(1, n ^ k);
        for i = 1:nTerms
                ps = unique(perms(indexSet(i, :)),'rows');
                for j = 1:mult(i) % add other permutations; mult(i) should be length(ps)
%                     Tv = Tv + kroneckerRight(v.', T(ps(j,:))); % Identity
                    Tv = Tv + kronMonomialSymmetrize(kroneckerRight(v.', T(ps(j,:))), n, k).'; % Can also symmetrize each thing
                end
        end
        
        Tv = Tv.';
        
    case 3 %% Now symmetrize and multiply by multiplicity rather than using perms
        % Get a list of indices
        [indexSet, mult] = findCombinations(m, k);
        nTerms = size(indexSet, 1);
        %
        n = size(T{1}, 2);
        
        Tv = zeros(1, n ^ k);
        for i = 1:nTerms
            Tv = Tv + mult(i) * kronMonomialSymmetrize(kroneckerRight(v.', T(indexSet(i, :))), n, k).'; % Can also symmetrize each thing
        end
        
        Tv = Tv.';
    
end
end

function indexSet = mkIndices(m, k)
%  Returns all combinations of m natural numbers that sum to k.
%
%  indexSet = mkIndices(m,k);
%
%  Used to compute calT_{m,k} in the function calTtimesv
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference:  Nonlinear balanced truncation model reduction for large-scale
%              polynomial systems, arXiv
%
%              See Algorithm 3.
%
%  Part of the NLbalancing repository.
%%
indexSet = [];

if m == 1
    indexSet = k;
    
else % fill the remaining indices through recursion
    
    for i = 1:k - m + 1
        J = mkIndices(m - 1, k - i);
        rJ = size(J, 1);
        
        indexSet = [indexSet; i * ones(rJ, 1) J];
    end
    
end

end

function [combinations, multiplicities] = findCombinations(m, k)
combinations = [];
multiplicities = [];
findCombinationsHelper(m, k, zeros(1, m), 1, 1);

    function findCombinationsHelper(m, k, combination, index, start)
        if index > m
            if k == 0
                combinations = [combinations; combination];
                multiplicities = [multiplicities; computeMultiplicity(combination)];
            end
        else
            for i = start:k
                combination(index) = i;
                findCombinationsHelper(m, k - i, combination, index + 1, i);
            end
        end
    end
end

% Calculate multiplicity correctly
function multiplicity = computeMultiplicity(comb)
n = length(comb);
multiplicity = factorial(n);
for i = unique(comb)
    count = sum(comb == i);
    multiplicity = multiplicity / factorial(count);
end
end

