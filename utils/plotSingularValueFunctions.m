function [] = plotSingularValueFunctions(sigma,c,zRange,n,maxDegree)
%  Plots polynomial approximations to singular value functions.
%
%   plotSingularValueFunctions(sigma,c,zRange,n,l)
%
%   Given an approximation to the singular value functions with the form
%
%       xi_i(z) = sigma(i) + sum_{k=1}^l diag(c{k})*z.^k
%
%   this function plots the above approximation for i=1,...,n at the points
%   specified in zRange.
%
%   The default value of n (the number of singular value functions plotted)
%   is determined from the length of sigma.
%
%   The default value of the maximum degree of the approximation, l, is
%   determined from the length of the cell array c.
%
%   Note that we explicitly assume the same range of z in all of the
%   coordinate directions.
%
%   Author: Jeff Borggaard, Virginia Tech
%           Modified by Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference: Nonlinear balanced trunction: Part 2--Model reduction on
%              manifolds, by Kramer, Gugercin, and Borggaard, arXiv.
%
%   Part of the NLbalancing repository.
%%
if (nargin<4)
    n = length(sigma);
end
if (nargin<5)
    maxDegree = length(c);
end

% Compute singular value function \xi; first term is diag(sigma).1,...
xi = sigma.*ones(n,length(zRange));% This ones() is not the 1 in the paper, it is just for plotting over zRange. sigma(i) does the job of Sigma . ones()
% ...then follows the sigma sum implemented as a for loop
for k=1:maxDegree
    xi = xi + c{k}.*repmat(zRange.^k,n,1);
end

figure
plot(zRange,xi);
legend(sprintfc('$\\xi_%d$',1:n), 'interpreter', 'latex')
xlabel('$z_{1,2}$', 'interpreter', 'latex');
ylabel('$\xi_{1,2}$', 'interpreter', 'latex');

if n == 2
    figure
    plot(zRange,xi);
    legend(sprintfc('$\\xi_%d$',1:n), 'interpreter', 'latex')
    xlabel('$z_{1,2}$', 'interpreter', 'latex');
    ylabel('$\xi_{1,2}$', 'interpreter', 'latex');
end



end