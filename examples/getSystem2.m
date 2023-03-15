function [A, B, C, N, g, f, h] = getSystem2()
%getSystem2  Generates a 2D polynomial system for testing energy functions.
%            System based on systems taken from [1,2].
%
%   Usage:  [A,B,C,N] = getSystem2()
%        or [~,~,~,~,~,f,g,h] = getSystem2()
%
%        \dot(x1) = -x1 + x2 - x2^2 + u1 + 2 x2 u1
%        \dot(x2) =     - x2        + u1
%              y1 =  x1
%
%   Excluding the higher order g terms gives the same system as in [1]. The
%   cell arrays f, g, and h can also be used as outputs.
%
%       References: [1] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki,
%                       “Nonlinear balanced truncation: part 1-computing energy
%                       functions,” Sep. 2022, doi: 10.48550/ARXIV.2209.07645.
%                   [2] Y. Kawano and J. M. A. Scherpen, “Model reduction by
%                       differential balancing based on nonlinear hankel operators,”
%                       IEEE Transactions on Automatic Control, vol. 62, no. 7,
%                       pp. 3293–3308, Jul. 2017, doi: 10.1109/tac.2016.2628201.
%%

A = [-1 1;
     0 -1];
N = [0 0 0 -1;
     0 0 0 0];
B = [1;
     1];
C = [1 1];
G1 = [0 2;
      0 0];
G2 = [0 0 0 0;
      0 0 0 0];

f = {A, N};
g = {B, G1, G2};
h = {C};

end
