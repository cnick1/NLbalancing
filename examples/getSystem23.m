function [f, g, h] = getSystem23(linear)
%getSystem23  Returns the 3D nonlinear model from [1,2].
%
%   Usage:  [f,g,h] = getSystem23()
%
%   Inputs:     linear - boolean, whether to use the linear model from [3]
%                        or nonlinear model from [1,2]
%
%   Outputs:     f,g,h - Cell arrays containing the polynomial coefficients
%                        for the drift, input, and output
%
%   Description: The dynamics correspond to the input-output system
%
%       ẋ₁ = −x₁ + 20 x₁ x₃ + u,
%       ẋ₂ = −2 x₂ + 20 x₂ x₃ + u,
%       ẋ₃ = −5 x₃ + u,
%        y = x₁ + x₂ + x₃,
%
%   If the option linear is enabled, the model from [3, Section 5.6.1] is
%   returned instead. This model replaces the nonlinear interaction between
%   x₃ and x₁ & x₂ with a linear interaction. The model from [3] inspired
%   the model in [1,2].
%
%       ẋ₁ = −x₁ + 100 x₃ + u,
%       ẋ₂ = −2 x₂ + 100 x₃ + u,
%       ẋ₃ = −5 x₃ + u,
%        y = x₁ + x₂ + x₃.
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


n = 3;
x = sym('x', [1, n]).'; syms(x);

if linear
    fsym = [-x1 + 100*x3;
        -2*x2 + 100*x3;
        -5*x3];
    gsym = [1;1;1];
    hsym = x1+x2+x3;
else
    fsym = [-x1 + 20*x1*x3;
        -2*x2 + 20*x2*x3;
        -5*x3];
    gsym = [1;1;1];
    hsym = x1+x2+x3;
end


[f, g, h] = approxPolynomialDynamics(fsym, gsym, hsym, x, 3);

end
