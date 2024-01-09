function [] = runExample9()
%runExample9 Runs the Allen-Cahn example.
%
%   Usage:  [v,w] = runExample9()
%
%   Inputs:
%
%   Outputs:
%
%   Background: Based on p34.m from [1].
%
%   Reference: [1] L. N. Trefethen, Spectral methods in MATLAB. Society
%              for Industrial and Applied Mathematics, 2000.
%              doi: 10.1137/1.9780898719598.
%
%   Part of the NLbalancing repository.
%%

fprintf('Running Example 9\n')

% Get system
% for y0 = [-.5 -.25 0 .25 .5]
eps = 0.01; N = 64; y0 = -0.5;
[f, g, ~, ~, y] = getSystem9(eps, N, y0);
fprintf("Maximum eigenvalue of A is %f; should be zero I think.\n",max(eigs(full(f{1}),N+1)))
% end

% Reference configuration (@ origin) -> v = v+vref
vref = tanh((y-y0)/sqrt(2*eps));

%% Solve PDE by Euler formula and plot results:
% Initial condition
v = .53*y + .47*sin(-1.5*pi*y);
% v = tanh((y-(0.25))/sqrt(2*eps));

% Shift initial condition 
v=v-vref;
% Time-stepping
dt = min([.0001,50*N^(-4)/eps]); t = 0; 
tmax = 100; tplot = 2; nplots = round(tmax/tplot);
plotgap = round(tplot/dt); dt = tplot/plotgap
xx = -1:.025:1; vv = polyval(polyfit(y,v+vref,N),xx);
plotdata = [vv; zeros(nplots,length(xx))]; tdata = t;
for i = 1:nplots
    fprintf('%i',i)
    for n = 1:plotgap
        % t = t+dt; v = v + dt*(f{1}*v + f{2}*kron(v,v) + f{3}*kron(v,kron(v,v)));    % Euler
        t = t+dt; v = v + dt*(f{1}*v + f{2}*kron(v,v) - v.^3);    % Euler
    end
    vv = polyval(polyfit(y,v+vref,N),xx);
    plotdata(i+1,:) = vv; tdata = [tdata; t];
end
clf, subplot('position',[.1 .4 .8 .5])
mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -1 1]),
view(-60,55), colormap(1e-6*[1 1 1]); xlabel x, ylabel t, zlabel u

%% Analytical solution
vref2 = tanh((xx+y0)/sqrt(2*eps));

hold on
plot3(xx,100*ones(size(xx)),vref2,"*")

return;

figure
hold on
plot(xx,plotdata(end,:))
plot(xx,vref2,"--")


figure
hold on
plot(y,v,"*")
plot(xx,vref2,"--")




function [D,x] = cheb(N)
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries
end

end
