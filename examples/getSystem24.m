function [f, g, h] = getSystem24(linear)
%getSystem24 Returns the 2D nonlinear model based on [1-3].
%
%   Usage:  [f,g,h] = getSystem24()
%
%   Inputs:     linear - boolean, whether to use the linear or nonlinear model
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The dynamics correspond to the input-output system
%
%       ẋ₁ = −2 x₁ + 20 x₁ x₂ + u,
%       ẋ₂ = −5 x₂ + u,
%        y = x₁ + x₂
%
%   If the option linear is enabled, the nonlinear interaction between x₃
%   and x₁ & x₂ is replaced with a linear interaction.
%
%       ẋ₁ = −2 x₁ + 100 x₂ + u,
%       ẋ₂ = −5 x₂ + u,
%        y = x₁ + x₂.
%
%   This is a 2D version of the models in getSystem23(), which come from [1-3].
%
%   References: [1] S. E. Otto, A. Padovan, and C. W. Rowley, "Optimizing
%                   oblique projections for nonlinear systems using
%                   trajectories," SIAM Journal on Scientific Computing,
%                   vol. 44, no. 3, pp. A1681–A1702, Jun. 2022, doi:
%                   10.1137/21m1425815.
%               [2] S. E. Otto, A. Padovan, and C. W. Rowley, "Model
%                   reduction for nonlinear systems by balanced truncation
%                   of state and gradient covariance,” SIAM Journal on
%                   Scientific Computing, vol. 45, no. 5, pp. A2325–A2355,
%                   Sep. 2023, doi: 10.1137/22m1513228.
%               [3] P. Holmes, J. L. Lumley, G. Berkooz, and C. W. Rowley,
%                   Turbulence, coherent structures, dynamical systems and
%                   symmetry. Cambridge University Press, 2012. doi:
%                   10.1017/cbo9780511919701.
%
%%

if nargin < 1
    linear = false;
end


n = 2;
x = sym('x', [1, n]).'; syms(x);

if linear
    fsym = [-2*x1 + 100*x2;
        -5*x2];
    gsym = [1;1];
    hsym = x1+x2;
else
    fsym = [-2*x1 + 20*x1*x2;
        -5*x2];
    gsym = [1;1];
    hsym = x1+x2;
end


[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, 3);

end
