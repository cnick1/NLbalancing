function [f, g, h, xdot] = getSystem22(degree)
%getSystem22 4D pendulum on a cart model.
%
%   Usage: [f, g, xdot] = getSystem22()
%
%   Background: The derivation is based on [1-3]. The result is a
%               control-affine nonlinear system, which we approximate with
%               Taylor expansions.
%
%   References: [1] https://www.mathworks.com/help/control/ref/lti.lqr.html
%               [2] https://www.mathworks.com/help/symbolic/derive-and-simulate-cart-pole-system.html.
%               [3] https://www.mathworks.com/help/mpc/ug/swing-up-control-of-a-pendulum-using-nonlinear-model-predictive-control.html
%
%   Part of the PPR repository.
%%
if nargin < 1 
    degree = 7;
end

%% Set up symbolic variables
syms t grav l I m M H V F theta(t) x(t) b u; X = sym('x', [1, 4]).';
assume([t grav l m M b] > 0); assume([I H V F],"real");

%% Set up equations of motion
eq1 = F - H - b*diff(x,1) == M*diff(x,2);                 % Sum of forces
eq2 = H == m*diff(x - l*sin(theta),2);                    % Sum of moments
mainEq1 = eliminate([eq1 eq2],H) == 0;

eq3 = -V*sin(theta) + m*grav*sin(theta) - H*cos(theta) == m*l*diff(theta,2) - m*diff(x,2)*cos(theta);
eq4 = V*l*sin(theta) + H*l*cos(theta) == I*diff(theta,t,2);
mainEq2 = eliminate([eq3 eq4],[H V]) == 0;

[ssEqs,states] = odeToVectorField([mainEq2 mainEq1]);
Ydot = matlabFunction(ssEqs,Vars=["t","Y","F","I","M","b","grav","l","m"]);

%% Specify physical quantities
Mval = 0.4801434687721852; mval = 0.18652319789448174; lval = 0.3574175621006708; 
bval = 0.1; Ival = 1/3 * mval * lval^3;  gval = 9;
xdot = @(X,u) Ydot(0,X,u,Ival,Mval,bval,gval,lval,mval);

%% Nonlinear equations of motion
vpa(xdot(X,u),2);

%% Linear equations of motion
dFdx = jacobian(xdot(X,u),X); vpa(subs(dFdx,X, [0;0;0;0]),2);
dFdu = jacobian(xdot(X,u),u); vpa(subs(dFdu,X, [0;0;0;0]),2);
 
%% Control-affine dynamics 
GofX = dFdu; FofX = simplify(xdot(X,u) - GofX*u);

%% Polynomial dynamics 
[f,g,~] = approxPolynomialDynamics(FofX,GofX,X,X,degree);
f{1} = full(f{1}); g{1} = full(g{1});
C = [1 0 0 0; 0 0 1 0]; h = {C};

% save() % Could write it to save to a .mat file so you only run once
end