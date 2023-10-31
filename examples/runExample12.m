function runExample12(degree)
%runExample12 Runs the example to test diagonalization. The system is a
%   polynomial approximation of the 2D model from Fujimoto and Scherpen
%   2001, 2005, 2010 [1-3].
%
%   Usage:  [v,w] = runExample12(degree)
%
%   Inputs:
%       degree - desired degree of the energy function approximation
%
%   Outputs:
%
%   The value of eta is set below.
%
%   References: [1] K. Fujimoto and J. M. A. Scherpen, “Model reduction
%                for nonlinear systems based on the differential
%                eigenstructure of Hankel operators,” in Proceedings of
%                the 40th IEEE Conference on Decision and Control (Cat.
%                No.01CH37228), IEEE, 2001. doi: 10.1109/cdc.2001.980322
%               [2] K. Fujimoto and J. M. A. Scherpen, “Nonlinear
%                input-normal realizations based on the differential
%                eigenstructure of Hankel operators,” IEEE Transactions
%                on Automatic Control, vol. 50, no. 1, pp. 2–18, Jan.
%                2005, doi: 10.1109/tac.2004.840476
%               [3] K. Fujimoto and J. M. A. Scherpen, “Balanced
%                realization and model order reduction for nonlinear
%                systems based on singular value analysis,” SIAM Journal
%                on Control and Optimization, vol. 48, no. 7, pp.
%                4591–4623, Jan. 2010, doi: 10.1137/070695332
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 12\n')
eta = 0;

if nargin < 1
    degree = 8;
end

[f, g, h] = getSystem12(degree - 1);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

fprintf("Beginning comparisons with Fujimoto/Scherpen 2001/2005/2010:\n") 

%% Comparison with Fujimoto 2001/2005 Example 1
x = sym('x', [1, 2]).'; syms(x);

LoSym =(36*x1^2 + 9*x2^2 + 18*x1^3*x2 + 18*x1*x2^3 + x1^6 + 6*x1^4*x2^2 + 9*x1^2*x2^4 + 4*x2^6)/ (1 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4);
[~,~,Lo] = approxPolynomialDynamics([1;1],[1;1],LoSym,x,degree);

fprintf("\n  - Comparing our energy function with Fujimoto/Scherpen 2001/2005 Example 1:\n") 
for i = 2:2:length(Lo)
    fprintf("    > Degree %i: the largest difference in w%i is %.1e;", ...
        i, i, norm(w{i} - kronMonomialSymmetrize(full(Lo{i}),2,i),'inf')) % Should be zero
    fprintf(" the largest entry in v%i is %.1e;\n", ...
        i, norm(v{i},'inf')) % Should be one for v2, zero else
end
fprintf("\n                             ->  Example 1 results match.\n\n")



%% Compute the input-normal transformation approximation
[sigma, Tin] = inputNormalTransformation(v, w, degree -1, true);

%% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
[sigmaSquared, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree - 1, true);

fprintf("\n  - Comparing our energy function with Fujimoto/Scherpen 2010 Example 5:\n") 

sigmaSquared;

s1Sym = 2*sqrt((9+x1^4)/(1+x1^4))*x1;
[~,~,s1] = approxPolynomialDynamics(1,1,s1Sym,x1,degree);
[~,~,s2] = approxPolynomialDynamics(1,1,s1Sym/2,x1,degree);
g1 = full(cell2mat(s1)); g1sq = conv(g1,g1);
g2 = full(cell2mat(s2)); g2sq = conv(g2,g2);

fprintf("    > Maximum distance between singular value functions squared:           %e \n",...
    norm([g1sq(1:degree-1); g2sq(1:degree-1)] - sigmaSquared,'inf'))
%% Plot the singular value functions

z = linspace(-1, 1, 51);
figure; hold on;

plot(z, sqrt(polyval(flip(sigmaSquared(1, :)), z)), 'LineWidth', 1.5); hold on
plot(z, polyval(flip(g1),z),'--', 'LineWidth', 1.5); hold on
plot(z, 2*sqrt((9+z.^4)./(1+z.^4)),'k--', 'LineWidth', 1.5)

plot(z, sqrt(polyval(flip(sigmaSquared(2, :)), z)), 'LineWidth', 1.5); hold on
plot(z, polyval(flip(g2),z),'--', 'LineWidth', 1.5); hold on
plot(z, sqrt((9+z.^4)./(1+z.^4)),'k--', 'LineWidth', 1.5)

fprintf("\n                             ->  Example 5 results match.\n\n")

end
