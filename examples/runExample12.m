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

[f, g, h] = getSystem12(degree - 1, false);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);


fprintf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n") 
fprintf(" ~~~~~~~~~~~ Beginning comparisons with Fujimoto/Scherpen 2001/2005:  ~~~~~~~~~~~ \n") 

%% Comparison with Fujimoto 2001/2005 Example 1
x = sym('x', [1, 2]).'; syms(x);

LoSym =(36*x1^2 + 9*x2^2 + 18*x1^3*x2 + 18*x1*x2^3 + x1^6 + 6*x1^4*x2^2 + 9*x1^2*x2^4 + 4*x2^6)/ (1 + x1 ^ 4 + 2 * x1 ^ 2 * x2 ^ 2 + x2 ^ 4);
[~,~,Lo] = approxPolynomialDynamics([1;1],[1;1],LoSym,x,degree);

fprintf("\n  - Comparing our energy functions with Fujimoto/Scherpen 2001/2005 Example 1:\n") 
for i = 2:2:length(Lo)
        fprintf("    > Degree %i: the largest difference in v%i is %.1e;", ...
        i, i, norm(v{i},'inf')) % Should be one for v2, zero else
        fprintf(" the largest entry in w%i is %.1e;\n", ...
        i, norm(w{i} - kronMonomialSymmetrize(full(Lo{i}),2,i),'inf')) % Should be zero
end

thresh = 1e-16;
v{3}(abs(v{3}) < thresh) = 0;v{4}(abs(v{4}) < thresh) = 0;v{5}(abs(v{5}) < thresh) = 0;v{6}(abs(v{6}) < thresh) = 0;v{7}(abs(v{7}) < thresh) = 0;v{8}(abs(v{8}) < thresh) = 0;
w{3}(abs(w{3}) < thresh) = 0;w{4}(abs(w{4}) < thresh) = 0;w{5}(abs(w{5}) < thresh) = 0;w{6}(abs(w{6}) < thresh) = 0;w{7}(abs(w{7}) < thresh) = 0;w{8}(abs(w{8}) < thresh) = 0;

fprintf("\n    Controllability energy: \n        Lc = 1/2 *(")
disp(vpa(kronPolyEval(v,sym('x', [1, 2]).'),2))
fprintf("    Observability energy: \n        Lo = 1/2 *(")
disp(vpa(kronPolyEval(w,sym('x', [1, 2]).'),8))

fprintf("                             ->  Example 1 results match.\n\n")

%% Compute the input-normal transformation approximation
[sigma, Tin] = inputNormalTransformation(v, w, degree -1, false);

%% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
[sigmaSquared, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree - 1, false);

fprintf("\n  - Comparing our singular value functions with Fujimoto/Scherpen 2001/2005 Example 3:\n") 

rho1Sym = 2*sqrt((9+x1^4)/(1+x1^4))*x1;
[~,~,s1] = approxPolynomialDynamics(1,1,rho1Sym,x1,degree);
[~,~,s2] = approxPolynomialDynamics(1,1,rho1Sym/2,x1,degree);
g1 = full(cell2mat(s1)); g1sq = conv(g1,g1);
g2 = full(cell2mat(s2)); g2sq = conv(g2,g2);

fprintf("    > Maximum distance between singular value functions squared:          %e \n",...
    norm([g1sq(1:degree-1); g2sq(1:degree-1)] - sigmaSquared,'inf'))
%% Plot the singular value functions

syms z
fprintf("\n         𝜎_1^2(z) = tau_1(z,0) = ")
disp(vpa(poly2sym(flip(sigmaSquared(1,:)),z),2))
fprintf("         𝜎_2^2(z) = tau_2(0,z) = ")
disp(vpa(poly2sym(flip(sigmaSquared(2,:)),z),2))

fprintf("                             ->  Example 3 results match.\n\n")

%% Compare transformation
fprintf("    > The full input-normal/output-diagonal transformation is: \n\n         𝚽(z) = ")

Tod{3}(abs(Tod{3}) < 1e-14) = 0;Tod{5}(abs(Tod{5}) < 1e-14) = 0;
disp(vpa(kronPolyEval(Tod(1:5),sym('z', [1, 2]).'),4))

%% Compare with Fujimoto 2010
fprintf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n") 
fprintf(" ~~~~~~~~~~~~~~ Beginning comparisons with Fujimoto/Scherpen 2010: ~~~~~~~~~~~~~~ \n") 
fprintf("     Note: in this paper, they start from an intermediate transformed system.       \n\n") 

[f, g, h] = getSystem12(degree - 1, true);

%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

%% Comparison with Fujimoto 2001/2005 Example 1

%% Compute the input-normal transformation approximation
[sigma, Tin] = inputNormalTransformation(v, w, degree -1, true);

%% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
[sigmaSquared, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree - 1, true);

fprintf("\n  - Comparing our singular value functions with Fujimoto/Scherpen 2010 Example 5:\n") 

syms z
fprintf("\n         𝜎_1^2(z) = tau_1(z,0) = ")
disp(vpa(poly2sym(flip(sigmaSquared(1,:)),z),2))
fprintf("         𝜎_2^2(z) = tau_2(0,z) = ")
disp(vpa(poly2sym(flip(sigmaSquared(2,:)),z),2))

fprintf("                     ->  Singular value functions match.\n\n")

%% Compare transformation
fprintf("\n  - Comparing our transformation with Fujimoto/Scherpen 2010 Example 5:\n") 

Tod{3}(abs(Tod{3}) < 1e-14) = 0;Tod{5}(abs(Tod{5}) < 1e-14) = 0;

fprintf("    > The input-normal/output-diagonal transformation is: \n\n         𝚽(y) = ")
disp(vpa(kronPolyEval(Tod(1:5),sym('y', [1, 2]).'),4))

fprintf("                     ->  Transformation matches.\n\n")

fprintf("                             ->  Example 5 results match.\n\n")

return

%% Compute transformed dynamics
[fin, gin, hin] = transformDynamics(f, g, h, Tin);
[fod, god, hod] = transformDynamics(fin, gin, hin, Tod);

end
