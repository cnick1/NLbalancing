function x = PhiBar(zbar,TinOd,sigmaSquared)
%PhiBar Return the balancing transformation x=Φbar(zbar) (at one point).
%
%   Usage:  x = PhiBar(zbar,TinOd,sigmaSquared)
%
%   Inputs:     zbar  - point at which to evaluate the transformation
%               TinOd - cell array containing input-normal/output-diagonal
%                       transformation coefficients
%        sigmaSquared - the coefficients of the squared singular value functions
%
%   Outputs:       x  - the value of x=Φbar(zbar)
%
%   Description: To Do
% 
%              x = Φ(z)
%                = T{1}*z + T{2}*(z⊗z) + ... + T{d}*(z...⊗z)
% 
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%

n = length(zbar); sigmabari = zeros(n,1);

sigmaSquared = flip(sigmaSquared,2); 

for i=1:n
    sigmabari(i) = real(polyval(sigmaSquared(i,:), zbar(i))^(1/2)); % real because the sigmasquared functions are sometimes negative (due to numerical error) so round to zero
end

% M = diag( sigmabari.^(-1/4) ); % negative exponent is like 1./x
% z = M*zbar; 

z = ( sigmabari.^(-1/2) ) .* zbar; % Use dot product to avoid forming n x n matrix

x = kronPolyEval(TinOd, z);

end