function [f1,f2] = plotDifferentialSingularValueFunctions(dataRange,varargin)
%  Plots polynomial approximations to singular value functions.
%
%   plotDifferentialSingularValueFunctions(sigma,c,zRange,n,l)
%
%   Author: Nick Corbin, UCSD
%
%   License: MIT
%
%   Reference:
%
%   Part of the NLbalancing repository.
%%
% [v, w] = runExample2(6,1,0,2,2,2,0,1);

% P = [5/4 3/4;
%     3/4 1/2];
% Q = [1/2 1/4;
%     1/4 1/4];

if nargin < 1
    dataRange = 1;
end

x1 = linspace(-dataRange, dataRange, 101);
x2 = x1;

[x1, x2] = meshgrid(x1, x2);

xPx = ((5/4+3*x2+2*x2.^2).*x1 + (3/4 + x2).*x2).*x1 + ((3/4 + x2).*x1 + 1/2*x2).*x2;
xPinvx = ((8).*x1 + (16*(-x2-3/4)).*x2).*x1 + ((16*(-x2-3/4)).*x1 + (16*(2*x2.^2+3*x2+5/4)).*x2).*x2;
xQx = ((1/2).*x1 + (1/4 - x2/3).*x2).*x1 + ((1/4 - x2/3).*x1 + (1/4 - 5*x2/9 + x2.^2/3).*x2).*x2;
% [8 , 16*(-x2-3/4) ;
%  16*(-x2-3/4) , 16*(2*x2.^2+3*x2+5/4)]

x=2;y=2;
((1/2).*x + (1/4 - y/3).*y).*x + ((1/4 - y/3).*x + (1/4 - 5*y/9 + y.^2/3)*y).*y;

% xPinvx = ((8).*x1 + (-12).*x2).*x1 + ((-12).*x1 + 20*x2).*x2;
% xQx = ((1/2).*x1 + (1/4).*x2).*x1 + ((1/4).*x1 + (1/4)*x2).*x2;

% figure
% contourf(x1,x2,.5*xPx)
% colorbar
if dataRange == 1
f1=figure('Name','Controllability/Past Energy Function')
contourf(x1,x2,.5*xPinvx,16,'w')
logMaxECTRB = log10(max(max(.5*xPinvx))); hold on;
contour(x1, x2, .5*xPinvx,[0,logspace(-3,ceil(logMaxECTRB),20)]./(10^(ceil(logMaxECTRB)-logMaxECTRB)))
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 16)
xticks([-dataRange 0 dataRange])
yticks([-dataRange 0 dataRange])
axis equal
load('utils\YlGnBuRescaled.mat')
colormap(flip(YlGnBuRescaled))
    caxis([0 80])


f2=figure('Name','Observability/Future Energy Function')
contourf(x1,x2,.5*xQx,16,'w')
logMaxEOBSV = log10(max(max(.5*xQx))); hold on;
contour(x1, x2, .5*xQx,[0,logspace(-3,ceil(logMaxEOBSV),20)]./(10^(ceil(logMaxEOBSV)-logMaxEOBSV)))
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
h=colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 16)
xticks([-dataRange 0 dataRange])
yticks([-dataRange 0 dataRange])
axis equal
load('utils\YlGnBuRescaled.mat')
colormap(flip(YlGnBuRescaled))
    caxis([0 1.5])
            set(h, 'ylim', [0 1.5])
else
   f1=figure('Name','Controllability/Past Energy Function')
contourf(x1,x2,.5*xPinvx,10,'w:','LineWidth',3)
% logMaxECTRB = log10(max(max(.5*xPinvx))); hold on;
% contour(x1, x2, .5*xPinvx,[0,logspace(-3,ceil(logMaxECTRB),20)]./(10^(ceil(logMaxECTRB)-logMaxECTRB)))
% xlabel('$x_1$', 'interpreter', 'latex');
% ylabel('$x_2$', 'interpreter', 'latex');
% colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
% set(gca, 'FontSize', 16)
xticks([])
yticks([])
axis equal
load('utils\YlGnBuRescaled.mat')
colormap(flip(YlGnBuRescaled))
    caxis([0 80])


f2=figure('Name','Observability/Future Energy Function')
contourf(x1,x2,.5*xQx,10,'w:','LineWidth',3)
% logMaxEOBSV = log10(max(max(.5*xQx))); hold on;
% contour(x1, x2, .5*xQx,[0,logspace(-3,ceil(logMaxEOBSV),20)]./(10^(ceil(logMaxEOBSV)-logMaxEOBSV)))
% xlabel('$x_1$', 'interpreter', 'latex');
% ylabel('$x_2$', 'interpreter', 'latex');
% colorbar('FontSize', 16, 'TickLabelInterpreter', 'latex')
% set(gca, 'FontSize', 16)
xticks([])
yticks([])
axis equal
load('utils\YlGnBuRescaled.mat')
colormap(flip(YlGnBuRescaled))
    caxis([0 1.5]) 
end

a = sqrt(81*x2.^4 + 396*x2.^3 + 844*x2.^2 + 900*x2+405);
lambda1 = 9/16 + 11*x2/18 + x2.^2/4 + a/36;
lambda2 = 9/16 + 11*x2/18 + x2.^2/4 - a/36;

% figure
% surf(x1,x2,sqrt(lambda1))
% colorbar
%
% figure
% surf(x1,x2,sqrt(lambda2))
% colorbar

%% z coordinates
z1 = -0.232*x1 - 0.143*x2 + 0.552e-3*x1.^2 + 0.682e-3*x1.*x2 + 0.110*x2.^2 - 0.523e-3*x1.*x2.^2 - 0.0214*x2.^3;
z2 =  0.214*x1 - 0.347*x2 + 0.552e-3*x1.^2 + 0.179e-2*x1.*x2 - 0.209*x2.^2 + 0.107e-2*x1.*x2.^2 - 0.273e-2*x2.^3;

% figure
% contourf(z1,z2,.5*xPinvx)
% colorbar
% xlim([-.2 .2])
% ylim([-.2 .2])
%
% figure
% contourf(z1,z2,.5*xQx)
% colorbar
% xlim([-.2 .2])
% ylim([-.2 .2])
%
% figure
% surf(z1,z2,sqrt(lambda1))
% shading interp
% colorbar
% xlim([-.2 .2])
% ylim([-.2 .2])
%
% figure
% surf(z1,z2,sqrt(lambda2))
% shading interp
% colorbar
% xlim([-.2 .2])
% ylim([-.2 .2])

a = sqrt(405);
sigma1 = sqrt(9/16 + a/36)
sigma2 = sqrt(9/16 - a/36)

end