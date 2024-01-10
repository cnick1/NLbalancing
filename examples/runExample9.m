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

%% Construct controller
% Get system
y0 = -.5; % Desired interface location

eps = 0.01; N = 32;
[f, B, ~, D, y] = getSystem9(eps, N, y0);
fprintf("Maximum eigenvalue of A is %f; should be zero I think.\n",max(eigs(full(f{1}),N+1)))

% Reference configuration (@ origin) -> v = v+vref
if isempty(y0)
    vref=zeros(N+1,1);
else
    vref = tanh((y-y0)/sqrt(2*eps));
end

B = B(:,linspace(1,N+1,5)); B(:,[1 5]) = []; m = size(B,2);
Q = 1; R = 1;


[ValueFun] = pqr(f, B, Q, R);
uOpenLoop = @(z) zeros(m,1);
uLinear = @(z) (- R * B.' * kronPolyDerivEval(ValueFun(1:2), z).' / 2);
uQuadratic = @(z) (- R * B.' * kronPolyDerivEval(ValueFun(1:3), z).' / 2);
uCubic = @(z) (- R * B.' * kronPolyDerivEval(ValueFun(1:4), z).' / 2);

controllers = {uOpenLoop, uLinear, uQuadratic, uCubic};

for idx = 1:4
    u = controllers{idx};

    %% Solve PDE by Euler formula and plot results:
    % Construct originial system dynamics
    D2 = D^2; D2([1 N+1],:) = zeros(2,N+1); % For boundary conditions

    % Initial condition
    v0 = .53*y + .47*sin(-1.5*pi*y);
    % v0 = tanh((y-(-0.125))/sqrt(2*eps*10));
    v = v0;

    % Time-stepping
    dt = min([.001,50*N^(-4)/eps]); t = 0;
    tmax = 100; tplot = 2; nplots = round(tmax/tplot);
    plotgap = round(tplot/dt); dt = tplot/plotgap;
    xx = -1:.025:1; vv = polyval(polyfit(y,v,N),xx);
    plotdata = [vv; zeros(nplots,length(xx))]; tdata = t;
    for i = 1:nplots
        fprintf('%i',i)
        for n = 1:plotgap
            t = t+dt; v = v + dt*(eps*D2*v + v - v.^3 + B*u(v-vref));    % Euler
        end
        vv = polyval(polyfit(y,v,N),xx);
        plotdata(i+1,:) = vv; tdata = [tdata; t];
    end
    figure, subplot('position',[.1 .4 .8 .5])
    mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -1 1]),
    view(-60,55), colormap(1e-6*[1 1 1]); xlabel x, ylabel t, zlabel u
    drawnow

end

return;

%% Solve PDE by Euler formula and plot results:
% Initial condition
v0 = .53*y + .47*sin(-1.5*pi*y);
% v0 = tanh((y-(-0.125))/sqrt(2*eps*10));
v = v0;

% Shift initial condition
v=v-vref;
% Time-stepping
dt = min([.001,50*N^(-4)/eps]); t = 0;
tmax = 100; tplot = 2; nplots = round(tmax/tplot);
plotgap = round(tplot/dt); dt = tplot/plotgap;
xx = -1:.025:1; vv = polyval(polyfit(y,v+vref,N),xx);
plotdata = [vv; zeros(nplots,length(xx))]; tdata = t;
for i = 1:nplots
    fprintf('%i',i)
    for n = 1:plotgap
        % t = t+dt; v = v + dt*(f{1}*v + f{2}*kron(v,v) + f{3}*kron(v,kron(v,v)));    % Euler
        t = t+dt; v = v + dt*(f{1}*v + f{2}*kron(v,v) - v.^3 + B*u(v));    % Euler
    end
    vv = polyval(polyfit(y,v+vref,N),xx);
    plotdata(i+1,:) = vv; tdata = [tdata; t];
end
figure, subplot('position',[.1 .4 .8 .5])
mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -1 1]),
view(-60,55), colormap(1e-6*[1 1 1]); xlabel x, ylabel t, zlabel u

%% Analytical solution

hold on
% plot3(xx,100*ones(size(xx)),polyval(polyfit(y,v0,N),xx),"ko")
plot3(xx,100*ones(size(xx)),polyval(polyfit(y,v+vref,N),xx),"*")

% load(fullfile('examples', 'systemData',sprintf('system9_equ ilibria_N=%i.mat',N)), 'd')
% plot3(xx,100*ones(size(xx)),polyval(polyfit(y,cell2mat(d(y0)),N),xx),"g*")

% view(180,0)

% figure
% hold on
% plot(xx,plotdata(end,:))
% plot(xx,polyval(polyfit(y,vref,N),xx),"--")
%
%
% figure
% hold on
% plot(y,v,"*")
% plot(xx,vref2,"--")




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
