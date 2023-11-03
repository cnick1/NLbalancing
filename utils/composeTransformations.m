function [Tout] = composeTransformations(T1,T2,varargin)
%composeTransformations Combines the transformations T1, T2, ... into one transformation 
%
%   Usage: [Tout] = composeTransformations(T1,T2)
%
%   Inputs:
%       T1,T2,... - cell arrays containing the polynomial transformation
%                   coefficients. At least two transformations are required,
%                   but more can be included and the function will apply
%                   the transformations recursively. 
%
%   Output:
%       Tout - cell arrays containing the polynomial coefficients for the
%              combined transformation.
%
%   TODO: figure out degree details
%
%   Background: Given a transformation T, compute the transformed dynamics.
%   TODO: Add more details here.
%
%   Authors: Nick Corbin, UCSD
%

ld = length(T1);
n = length(T1{1});

Tout = cell(size(T1));

for k = 1:ld
    Tout{k} = zeros(n,n^k);
%     for i=1:k % row by row
%         for j = 1:n
%             Tout{k}(j,:) = Tout{k}(j,:) + calTTv(T2, i, k, T1{i}(j,:).').';
%         end
%     end
        for i=1:k % should also work like this, just required calTTv and symmetrization functions to appropriately handle matrices
            Tout{k}= Tout{k} + calTTv(T2, i, k, T1{i}.').';
        end
end

if nargin > 2 % Compute the remaining ones by recursion
    [Tout] = composeTransformations(Tout,varargin{3},varargin{4:end});
end

end
