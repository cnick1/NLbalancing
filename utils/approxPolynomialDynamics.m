function [F,G,H] = approxPolynomialDynamics(f,g,h,x,degree)
%approxPolynomialDynamics Returns polynomial approximations to symbolic dynamics
%   
%   Usage: [F,G,H] = approxPolynomialDynamics(f,g,h,x,degree)
%
%   Inputs:
%       f,g,h   - symbolic expressions for the drift, input, and output of
%                 a dynamical system.
%       x       - vector of symbolic variables used in f,g,h; if not included,
%                 the default is "x = sym('x',[1,n])" where n is the length of f.
%       degree  - desired degree of the polynomial approximation to the
%                 dynamics.
%
%   Outputs:
%       F,G,H   - cell arrays containing the polynomial coefficients for
%                 the drift, input, and output. Computed using Yaylor
%                 approximations.
% 
%  Author: Nick Corbin, UCSD
%
%  License: MIT
%
%  Part of the NLbalancing repository.
%%

n = size(f,1);
m = size(g,2);
p = size(h,1);

if nargin < 5
    degree = 5;
    if nargin < 4
        x = sym('x',[1,n]);
        syms(x);
    end
end

% Initialize function outputs
F = cell(1,degree); G = cell(1,degree); H = cell(1,degree);

% Compute first entry of G
G{1} = double(taylor(g,x,'Order',1));

% Initialize previous Taylor expansion
prevTaylorExpansionF = 0;
prevTaylorExpansionG = taylor(g,x,'Order',1);
prevTaylorExpansionH = 0;

xn = 1;

% Compute each degree of taylor approximation, and iterate through grabbing
% the coefficients and putting them in a matrix
for d=1:degree
    F{d} = sparse(n,n^d);
    G{d+1} = sparse(n,m*n^(d));
    H{d} = sparse(p,n^d);
    
    currentTaylorExpansionF = taylor(f,x,'Order',d+1);
    currentTaylorExpansionG = taylor(g,x,'Order',d+1);
    currentTaylorExpansionH = taylor(h,x,'Order',d+1);
    
    taylorExpansionF = currentTaylorExpansionF - prevTaylorExpansionF;
    taylorExpansionG = currentTaylorExpansionG - prevTaylorExpansionG;
    taylorExpansionH = currentTaylorExpansionH - prevTaylorExpansionH;
    
    xn = kron(xn, x);
    
    % Construct F cell array
    for i=1:n % i is ith equation of n
        [Cf,Tf] = coeffs(taylorExpansionF(i),x);
        
        % Map T to xn
        locsF = map_symbolic_vectors(Tf, xn);
        
        %  Apply that mapping to C, write to coefficient
        try
            F{d}(i,locsF) = double(Cf);
        catch
            warning("Dynamics contain an affine contribution; discarding it")
            F{d}(i,locsF) = double(Cf(1:end-1));
        end
    end
    
    % Construct G cell array
    for i=1:n % i is ith equation of n
        for j=1:m
            [Cg,Tg] = coeffs(taylorExpansionG(i,j),x);
            
            % Map T to xn
            locsG = map_symbolic_vectors(Tg, xn);
            
            % Map to kron(x,u)
            locsG = ((locsG - 1)*m + 1) + j - 1;
            
            %  Apply that mapping to C, write to coefficient
            G{d+1}(i,locsG) = double(Cg);
        end
    end
    
    % Construct H cell array
    for i=1:p % i is ith equation of p
        [Ch,Th] = coeffs(taylorExpansionH(i),x);
        
        % Map T to xn
        locsH = map_symbolic_vectors(Th, xn);
        
        %  Apply that mapping to C, write to coefficient
        H{d}(i,locsH) = double(Ch);
    end
    
    prevTaylorExpansionF = currentTaylorExpansionF;
    prevTaylorExpansionG = currentTaylorExpansionG;
    prevTaylorExpansionH = currentTaylorExpansionH;
end

end


function mapping = map_symbolic_vectors(T, xn)
% Initialize an empty cell array to store the mapping
mapping = zeros(length(T),1);

% Loop through the elements of T and populate the mapping
for i = 1:length(T)
    [~,loc] = ismember(T(i), xn);
    mapping(i) = loc;
end
mapping(mapping == 0) = [];
end

