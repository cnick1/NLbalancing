function [A,B,C,N,Q] = getSystem7()
%getSystem7  Generates a quadratic bilinear system for testing energy functions.
%            System taken from [1]
%
%        \dot(x1) = -x1 + x2 - x2^2 + (1 + 2 x2) u1
%        \dot(x2) =     - x2        + u1
%              y1 =  x1 
%       References:
%%

  A = [-1  1;
        0 -1];
  N = [0 0 0 -1;
       0 0 0 0];
  B = [1;
       1];
  C = [1 0];
  Q = [0 2; 
      0 0];

end