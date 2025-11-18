function runExample17_checkDiag(n)
%runExample17 Runs the example to test diagonalization.
%
%   Usage:  [v,w] = runExample17_checkDiag()
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%
fprintf('Running Example 17\n')
close all;

% Create a vec function for readability
vec = @(X) X(:);

if nargin < 1
    n = 16;
end

eta = 0;

degree = 4;
[f, g, h] = getSystem17(degree - 1, n / 2);

% Compute the energy functions
[v] = approxPastEnergy(f, g, h, eta=eta, degree=degree, verbose=true);
[w] = approxFutureEnergy(f, g, h, eta=eta, degree=degree, verbose=true);

% Compute the input-normal/output-diagonal transformation approximation, also giving the squared singular value functions
[sigmaSquared, TinOd, vbar, wbar] = inputNormalOutputDiagonalTransformation(v, w, degree=degree-1, verbose=true);

%% Compute norms to evaluate diagonalization performance
% Save to file
fprintf('Writing data to plots/example17_norms_n%i.dat \n',n);
% fileID = 1;
fileID = fopen(sprintf('plots/example17_norms_n%i.dat',n), 'w');
fprintf(fileID, '# Table II Data \n');
fprintf(fileID, '# Output-Diagonalization Metrics \n');
%print the header
fprintf(fileID, '                                quantity                           &  2-norm  & infty-norm \n');
quantity = vbar{2}-vec(speye(n));
fprintf(fileID, '$\\left\\lVert \\tilde\\bv_2 - \\text{vec}(\\bI_n)\\right\\rVert$                      &   %1.1e    &   %1.1e \n', norm(quantity,2), norm(quantity,'inf'));
quantity = vbar{4};
fprintf(fileID, '$\\left\\lVert \\tilde\\bv_4\\right\\rVert$                              &   %1.1e    &   %1.1e \n', norm(quantity,2), norm(quantity,'inf'));
quantity = wbar{2}; quantity(linspace(1,n^2,n)) = 0; % make quantity the off-diagonal elements of w2
fprintf(fileID, '$\\left\\lVert \\tilde\\bw_2 - \\text{diag}(\\tilde\\bw_2)\\right\\rVert$   &   %1.1e    &   %1.1e \n', norm(quantity,2), norm(quantity,'inf'));
quantity = wbar{4}; quantity(linspace(1,n^4,n)) = 0; % make quantity the off-diagonal elements of w4
fprintf(fileID, '$\\left\\lVert \\tilde\\bw_4 - \\text{diag}(\\tilde\\bw_4)\\right\\rVert$   &   %1.1e    &   %1.1e \n', norm(quantity,2), norm(quantity,'inf'));


%% Shrink v4,w4 to 3 dimensions
v4=vec(vecnorm(reshape(v{4},n,n,n,n),2,4));
w4=vec(vecnorm(reshape(w{4},n,n,n,n),2,4));

v4 = v4/max(v4);
w4 = w4/max(w4);

% Save to file
fprintf('Writing data to plots/example17_projected3D-OGv4w4_n%i.dat \n',n);
fileID = fopen(sprintf('plots/example17_projected3D-OGv4w4_n%i.dat',n), 'w');
fprintf(fileID, '# Table I Data\n');
fprintf(fileID, '# OG v4 and w4 projected down to 3 dimensions; n=%i \n',n);
%print the header
fprintf(fileID, 'v4              & w4 \n');
for i = 1:length(v4)
    fprintf(fileID, '%12.6e    & %12.6e \n', v4(i), w4(i));
end
fclose(fileID);

vb4=vec(vecnorm(reshape(vbar{4},n,n,n,n),2,4));
wb4=vec(vecnorm(reshape(wbar{4},n,n,n,n),2,4));

vb4 = vb4/max(wb4);
wb4 = wb4/max(wb4);

% Save to file
fprintf('Writing data to plots/example17_projected3D-v4w4_n%i.dat \n',n)
fileID = fopen(sprintf('plots/example17_projected3D-v4w4_n%i.dat',n), 'w');
fprintf(fileID, '# Table I Data\n');
fprintf(fileID, '# v4 and w4 projected down to 3 dimensions; n=%i \n',n);
%print the header
fprintf(fileID, 'v4              & w4 \n');
for i = 1:length(v4)
    fprintf(fileID, '%12.6e    & %12.6e \n', vb4(i), wb4(i));
end
fclose(fileID);

%% Shrink v4,w4 to 2 dimensions
v4=vec(vecnorm(reshape(v{4},n,n,n,n),2,4));
w4=vec(vecnorm(reshape(w{4},n,n,n,n),2,4));

v4=vec(vecnorm(reshape(v4,n,n,n),2,3));
w4=vec(vecnorm(reshape(w4,n,n,n),2,3));

v4 = v4/max(v4);
w4 = w4/max(w4);

% Save to file
fprintf('Writing data to plots/example17_projected2D-OGv4w4_n%i.dat \n',n);
fileID = fopen(sprintf('plots/example17_projected2D-OGv4w4_n%i.dat',n), 'w');
fprintf(fileID, '# Table I Data\n');
fprintf(fileID, '# OG v4 and w4 projected down to 2 dimensions; n=%i \n',n);
%print the header
fprintf(fileID, 'v4              & w4 \n');
for i = 1:length(v4)
    fprintf(fileID, '%12.6e    & %12.6e \n', v4(i), w4(i));
end
fclose(fileID);

vb4=vec(vecnorm(reshape(vbar{4},n,n,n,n),2,4));
wb4=vec(vecnorm(reshape(wbar{4},n,n,n,n),2,4));

vb4=vec(vecnorm(reshape(vb4,n,n,n),2,3));
wb4=vec(vecnorm(reshape(wb4,n,n,n),2,3));

vb4 = vb4/max(wb4);
wb4 = wb4/max(wb4);

% Save to file
fprintf('Writing data to plots/example17_projected2D-v4w4_n%i.dat \n',n)
fileID = fopen(sprintf('plots/example17_projected2D-v4w4_n%i.dat',n), 'w');
fprintf(fileID, '# Table I Data\n');
fprintf(fileID, '# v4 and w4 projected down to 2 dimensions; n=%i \n',n);
%print the header
fprintf(fileID, 'v4              & w4 \n');
for i = 1:length(vb4)
    fprintf(fileID, '%12.6e    & %12.6e \n', vb4(i), wb4(i));
end
fclose(fileID);

%% Some plots

for j=2:2:length(vbar)
    figure
    heatmap(projectTensor(reshape(vbar{j},n*ones(1,j)), 2, 2))
end


for j=2:2:length(wbar)
    figure
    heatmap(projectTensor(reshape(wbar{j},n*ones(1,j)), 2, 2))
    
    figure
    heatmap(projectTensor2(reshape(wbar{j},n,n,[]), 2))
end

figure
wbar{j}(linspace(1,n^j,n)) = 0;
heatmap(projectTensor(reshape(wbar{j},n*ones(1,j)), 2, 2))
end

function A = projectTensor(A, p, rk)
% Iterate while contracting along fiber direction 1 until desired tensor rank
while length(size(A)) > rk
    A = squeeze(vecnorm(A,p,1));
end
end


function A = projectTensor2(A, p)
% Reshape A as an n by n by [] tensor; then just norm along the third direction
A = squeeze(vecnorm(A,p,3));
end

