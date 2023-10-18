function [A, B, C, F2, f, g, h] = getSystem11(nFterms, m, L, varargin)
%getSystem11  Generates a 2D inverted pendulum model with degree ell 
%             polynomial sin approximation, up to degree 9.
%   
%   Usage:  [A,B,C,N] = getSystem11(ell)
%        or [~,~,~,~,~,f,g,h] = getSystem11(ell)
%
%              \dot{x}_1 & = x_2                                                                                                                     \\
%              \dot{x}_2 & = \frac{3T}{m\ell^2} + \frac{3g}{2\ell} (x_1-\frac{x_1^3}{6}+\frac{x_1^5}{120}-\frac{x_1^7}{5040} + \frac{x_1^9}{362880})
%                   y    & = x_1
%
%       References: [1] 
%%

if nargin < 3
    if nargin < 2
        if nargin < 1
            nFterms = 7;
        end
        m = 1;
    end
    L = 1;
end

gravity = 9.81;

A = [0, 1; 3*gravity/(2*L), 0];
F2 = sparse(2,2^2); 
% F3 = sparse(2,2^3); F3(2,1) = -gravity/(4*L);
% F4 = sparse(2,2^4); 
% F5 = sparse(2,2^5); F5(2,1) = gravity/(80*L);
% F6 = sparse(2,2^6); 
% F7 = sparse(2,2^7); F7(2,1) = -gravity/(3360*L);
% F8 = sparse(2,2^8); 
% F9 = sparse(2,2^9); F9(2,1) = gravity/(241920*L);

f = {A, F2};
for i = 1:(nFterms-1)/2
    f{end + 1} = sparse(2,2^(2*i+1));
    f{end}(2,1) = 3*gravity/(2*L) * (-1)^i/factorial(2*i+1);
    f{end + 1} = sparse(2,2^(2*i+2));
end
f = f(1:nFterms);

B = [0; 3/(m*L^2)]; % Torque control at pivot point

C = [1, 0]; % Measure angle of pendulum


% f = {A, F2, F3, F4, F5, F6, F7, F8, F9}; f = f(1:nFterms);
g = {B};
h = {C};

end
