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
y0 = .5; % Desired interface location

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
Q2 = .1; Q3 = sparse((N+1)^3,1) ; Q4 = sparse(linspace(1,(N+1)^4,N+1),1,4);
q = {[],Q2,Q3,Q4};
R = 1;


[ValueFun] = pqr(f, B, q, R);
uOpenLoop = @(z) zeros(m,1);
uLinear = @(z) (- R * B.' * kronPolyDerivEval(ValueFun(1:2), z).' / 2);
uCubic = @(z) (- R * B.' * kronPolyDerivEval(ValueFun(1:4), z).' / 2);

controllers = {uOpenLoop, uLinear, uCubic};

for idx = 1:3
    u = controllers{idx};
    Lagrangian = zeros(100001,1);

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
            xbar = v-vref;
            Ux = u(xbar);
            Lagrangian(n+(i-1)*plotgap) = 1/2*(xbar.'*Q2*xbar + Ux.'*R*Ux + 4*sum(xbar.^4)); % hardcoded v.^4 instead of Q4 for speed
            t = t+dt; v = v + dt*(eps*D2*v + v - v.^3 + B*Ux);    % Euler
        end
        vv = polyval(polyfit(y,v,N),xx);
        plotdata(i+1,:) = vv; tdata = [tdata; t];
    end
    figure, subplot('position',[.1 .4 .8 .5])
    mesh(xx,tdata,plotdata), grid on, axis([-1 1 0 tmax -1 1]),
    view(-60,55), colormap([0 0 0]); xlabel x, ylabel t, zlabel u
    drawnow

    % Compute performance Index (cost)
    performanceIndex = trapz((0:dt:tmax), Lagrangian);
    fprintf("\n\n   The performance index is %f\n\n",performanceIndex)

end




exportgraphics(gcf,'plots/example9_cubic.pdf', 'ContentType', 'vector')
close
exportgraphics(gcf,'plots/example9_linear.pdf', 'ContentType', 'vector')
close
exportgraphics(gcf,'plots/example9_openloop.pdf', 'ContentType', 'vector')
close



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
