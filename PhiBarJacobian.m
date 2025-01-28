function J = PhiBarJacobian(zbar,TinOd,sigmaSquared)
%PhiBarJacobian Return the balancing transformation x=Φbar(zbar) (at one point).
%
%   Usage:  J = PhiBarJacobian(zbar,TinOd,sigmaSquared)
%
%   Inputs:     zbar  - point at which to evaluate the transformation
%               TinOd - cell array containing input-normal/output-diagonal
%                       transformation coefficients
%        sigmaSquared - the coefficients of the squared singular value functions
%
%   Outputs:       J  - the value of the Jacobian of the balancing transformation
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

n = length(zbar); sigmabari = zeros(n,1); dsigmabari = zeros(n,1);

sigmaSquared = flip(sigmaSquared,2); 
dsigmaSquared = zeros(size(sigmaSquared));

for i=1:n
    sigmabari(i) = real(polyval(sigmaSquared(i,:), zbar(i))^(1/2)); % real because the sigmasquared functions are sometimes negative (due to numerical error) so round to zero
    
    dsigmaSquared(i,2:end) = polyder(sigmaSquared(i,:));
    dsigmabari(i) = polyval(dsigmaSquared(i,:), zbar(i)) / (2*sigmabari(i));
end

dphi_ii = ( sigmabari.^(1/2) - 1/2 * zbar .* sigmabari.^(-1/2) .* dsigmabari ) ./ sigmabari; 

% J = jcbn(TinOd, dphi .* zbar);
J = jcbn(TinOd, sigmabari.^(-1/2).*zbar)  * diag(dphi_ii);


end