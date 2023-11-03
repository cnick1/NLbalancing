function [Tout] = composeTransformations(T1,T2,varargin)
%composeTransformations Transforms the dynamics f, g, h by T.
%
%   Usage: [ft, gt, ht] = composeTransformations(f, g, h, T)
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

ld = length(T1);
n = length(T1{1});

Tout = cell(size(T1));

for k = 1:ld
    Tout{k} = zeros(n,n^k);
%     for i=1:k
%         for j = 1:n
%             Tout{k}(j,:) = Tout{k}(j,:) + calTTv(T2, i, k, T1{i}(j,:).').';
%         end
%     end
        for i=1:k
            Tout{k}= Tout{k} + calTTv(T2, i, k, T1{i}.').';
        end
end

if nargin > 2 % Compute the remaining ones by recursion
    [Tout] = composeTransformations(Tout,varargin{3},varargin{4:end});
end

end
