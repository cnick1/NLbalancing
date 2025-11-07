function Tinv = transformationInverse(T)
%transformationInverse Return the coefficients of the inverse of the transformation.
%
%   Usage:  Tinv = transformationInverse(T)
%
%   Inputs: T - cell array containing transformation coefficients
%
%   Description: We approximate the transformation
%   as a polynomial, so
%
%              x = Φ(z)
%                = T{1}*z + T{2}*(z⊗z) + ... + T{d}*(z...⊗z)
%
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
Tinv = cell(size(T));

n = size(T{1}, 1); In = speye(n);
Tinv{1} = inv(T{1}); % Should be able to do explicitly

for i=2:length(T)
    Tinv{i} = zeros(size(T{i}));
    for j = 1:(i-1) % index for P_i
        for idx = 1:n
            Tinv{i}(idx,:) = Tinv{i}(idx,:) - calTTv(T,j,i,Tinv{j}(idx,:).').';
        end
    end
    % Now invert on right
    for idx = 1:n
        Tinv{i}(idx,:) = calTTv(Tinv,i,i,Tinv{i}(idx,:).').';
    end
end

end