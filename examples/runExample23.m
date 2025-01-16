function runExample23(degree,linear)
%runExample23 Runs the 3D model from Clancy Rowley's talk
%
%   Usage:  runExample23(degree)
%
%   Inputs:
%       degree - desired degree of the energy function approximation
%       linear - boolean, whether to use the linear model or nonlinear
%                model from the talk
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
set(groot,'defaultLineLineWidth',1.5)

fprintf('Running Example 23\n')

if nargin < 2
    if nargin < 1
        degree = 10;
    end
    linear = false;
end

[f, g, h] = getSystem23(linear);

%%  Compute the energy functions
fprintf("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~ Computing energy functions:  ~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
fprintf(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

eta = 0;
[v] = approxPastEnergy(f, g, h, eta, degree, false);
[w] = approxFutureEnergy(f, g, h, eta, degree, false);

% Print some results
n = 3;
x = sym('x', [1, n]).'; syms(x);

fprintf("\n    Controllability energy: \n        Lc = 1/2 *(")
disp(vpa(kronPolyEval(v, x), 3))
fprintf("    Observability energy: \n        Lo = 1/2 *(")
disp(vpa(kronPolyEval(w, x), 3))

%% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
tic
[sigmaSquared, Tod] = inputNormalOutputDiagonalTransformation(v, w, degree - 1, true);
fprintf("Input-normal/output-diagonal transformation took %f seconds. \n\n", toc)

%% Plot the squared singular value functions
fprintf("    Singular value functions (squared): \n\n")

% figure; hold on;
syms z
for i=1:n
    fprintf("         𝜎_%i^2(z) = tau_%i(z,0) = ",i,i)
    disp(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 4))
    fplot(vpa(poly2sym(flip(sigmaSquared(i, :)), z), 4),[-1 1])
end

set(gca, 'YScale', 'log')

return
%% Compute transformed dynamics
[ft, gt, ht] = transformDynamics(f, g, h, Tod);

% Compare original dynamics with transformed dynamics 
z = sym('z', [1, n]).'; syms(z); syms u;
F = @(x) kronPolyEval(f, x); G = @(x) (g{1} + kronPolyEval(g(2:end), x));
Ft = @(x) kronPolyEval(ft, x); Gt = @(x) (gt{1} + kronPolyEval(gt(2:end), x));

fprintf("\n    Original dynamics: xdot = \n")
disp(vpa(F(x) + G(x)*u, 3))
fprintf("\n    Transformed dynamics: zdot = \n")
disp(vpa(Ft(z) + Gt(z)*u, 3))


%% Simulate original and transformed systems; same input should give same output
x0 = [1 1 1].'/30; 

% [t1, X1] = ode45(@(t, x) F(x), [0, 5], x0*.5);
[t2, X2] = ode45(@(t, x) F(x), [0, 5], x0);

% figure
% subplot(2,1,1)
% plot(t1,X1)
% subplot(2,1,2)
% plot(t1,sum(X1,2))
% 
% figure
% subplot(2,1,1)
% plot(t2,X2)
% subplot(2,1,2)
% plot(t2,sum(X2,2))

figure
hold on
subplot(2,1,1); hold on;
plot(t2,X2)

% plot(t1,sum(X1,2))
subplot(2,1,2); hold on;
plot(t2,sum(X2,2))

z0 = Tod{1}\x0; z0.';
for i=1:10
    z0 = Tod{1}\(x0-kronPolyEval(Tod, z0)+Tod{1}*z0);
    z0.';
    kronPolyEval(Tod, z0);
end
kronPolyEval(Tod, z0)


% [t12, Z1] = ode45(@(t, x) Ft(x), [0, 5], z0*.5);
[t22, Z2] = ode45(@(t, x) Ft(x), [0, 5], z0);

X22 = zeros(size(Z2));
for i=1:length(t22)
    X22(i,:) = kronPolyEval(Tod, Z2(i,:).');
end

subplot(2,1,1)
plot(t22,X22,'--')
subplot(2,1,2)
plot(t22,sum(X22,2),'--')

%% Truncate last state to obtain a reduced-order model
% for i=1:length(ft)
%     fr{i} = ft{i}(1:2,:);
% end
% 
% for i=1:length(gt)
%     fr{i} = ft{i}(1:2,:);
% end
% 
% Fr = @(x) kronPolyEval(fr, x); Gr = @(x) (gr{1} + kronPolyEval(gr(2:end), x));



end
