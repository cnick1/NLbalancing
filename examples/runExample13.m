function runExample13()
%runExample13 Runs the example to test diagonalization. The system is a
%   polynomial approximation of the 2D model from Gray and Scherpen [1].
%
%   Usage:  [v,w] = runExample13()
%
%   References: [1] W. S. Gray and J. M. A. Scherpen, “On the nonuniqueness
%               of singular value functions and balanced nonlinear
%               realizations,” Systems & Control Letters, vol. 44, no. 3,
%               pp. 219–232, Oct. 2001, doi: 10.1016/s0167-6911(01)00144-x
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 13\n')
eta = 0;

[f, g, h] = getSystem13();

degree = 4;
%  Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta, degree, true);
[w] = approxFutureEnergy(f, g, h, eta, degree, true);

fprintf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
fprintf(" ~~~~~~~~~~~~~~~~ Beginning comparisons with Gray/Scherpen 2001: ~~~~~~~~~~~~~~~~ \n")

%% Comparison with Gray 2001 Example 1
fprintf("\n  - Comparing our energy function with Gray/Scherpen 2001 Example 2.1:\n")

v{3}(abs(v{3}) < 1e-14) = 0; v{4}(abs(v{4}) < 1e-14) = 0; w{3}(abs(w{3}) < 1e-14) = 0; w{4}(abs(w{4}) < 1e-14) = 0;

fprintf("\n    > Controllability energy: \n        Lc = 1/2 *(")
disp(vpa(kronPolyEval(v, sym('x', [1, 2]).'), 2))
fprintf("    > Observability energy: \n        Lo = 1/4 *(")
disp(vpa(kronPolyEval(w, sym('x', [1, 2]).') * 2, 2))

fprintf("                             ->  Energy functions match.\n\n")

%% Compute the input-normal transformation approximation
[sigma, Tin] = inputNormalTransformation(v, w, degree -1, true);

%% Compute the output-diagonal transformation approximation, also giving the squared singular value functions
[sigmaSquared, Tod] = outputDiagonalTransformation(v, w, Tin, diag(sigma), degree - 1, true);

fprintf("\n  - Comparing our singular value functions with Gray/Scherpen 2001 Example 2.1:\n")

sigmaSquared(sigmaSquared < 1e-14) = 0;

syms z
fprintf("\n         𝜎_1^2(z) = tau_1(z,0) = ")
disp(vpa(poly2sym(flip(sigmaSquared(1, :)), z), 2))
fprintf("         𝜎_2^2(z) = tau_2(0,z) = ")
disp(vpa(poly2sym(flip(sigmaSquared(2, :)), z), 2))

fprintf("\n                             ->  Singular value functions match up to 1e-14.\n\n")

%% Compare transformation
fprintf("\n  - Comparing our transformation with Gray/Scherpen 2001 Example 2.1:\n")

fprintf("    > The full input-normal/output-diagonal transformation is: \n\n         𝚽(z) = ")
disp(vpa(kronPolyEval(Tod, sym('z', [1, 2]).'), 2))

%% Gray transformation
z = sym('z', [1, 2]).'; syms(z);
% Tsym = [(-1 + sqrt(1+4*z1))/2; z2];
% [Tod2,~,~] = approxPolynomialDynamics(Tsym,[1;1],z1,z,3);
Tod2 = {[1 0; 0 1], [0 0 0 -1; 0 0 0 0], zeros(2, 2 ^ 3)};

fprintf("    > The transformation given in Gray/Scherpen 2001 is:     𝚽(z) = ")
% fprintf('%s \n', char(Tsym))
% fprintf("     which can be approximated via Taylor series as: \n\n         𝚽(z) = ")
fprintf('%s \n', char(vpa(kronPolyEval(Tod2, sym('z', [1, 2]).'), 2)))


[vtilde, wtilde] = transformEnergyFunctions(v, w, Tod2);
[vtilde, wtilde] = transformEnergyFunctions(vtilde, wtilde, {[1,-1; 1,1]./sqrt(2), zeros(2, 2 ^ 2), zeros(2, 2 ^ 3)});

thresh = 2e-14;
vtilde{3}(abs(vtilde{3}) < thresh) = 0; vtilde{4}(abs(vtilde{4}) < thresh) = 0; wtilde{3}(abs(wtilde{3}) < thresh) = 0; wtilde{4}(abs(wtilde{4}) < thresh) = 0;

fprintf("\n    > Controllability energy: \n        Lc = 1/2 *(")
disp(vpa(kronPolyEval(vtilde, sym('x', [1, 2]).'), 2))
fprintf("    > Observability energy: \n        Lo = 1/2 *(")
disp(vpa(kronPolyEval(wtilde, sym('x', [1, 2]).'), 2))

end
