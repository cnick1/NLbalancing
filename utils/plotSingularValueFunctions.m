function plotSingularValueFunctions(sigmaSquared)
%plotSingularValueFunctions Produce singular value function plot
%
%   Usage: plotSingularValueFunctions(sigmaSquared)
%
%   Inputs:
%       sigmaSquared - polynomial coefficient matrix for the squared
%                      singular value functons, which is produced as an 
%                      output of inputNormalOutputDiagonalTransformation() 
%                      and getBalancedThenReduceRealization().
%   
%   Description: The singular value functions provide the ranking of the
%   relative importance of state components, which can be used to decide
%   upon the truncation order for the reduced-order model. This function
%   plots the singular value functions on a local domain to show the user
%   this information.
%

n = size(sigmaSquared,1);
z = linspace(- 1, 1, 100001);
figure;
subplot(2,1,1); 
hold on; title("Singular Value Functions")
for i = 1:n
    temp = real(sqrt(polyval(flip(sigmaSquared(i, :)), z)));
    plot(z(temp>0), temp(temp>0),'DisplayName',sprintf('\\sigma_{%i}',i),'LineWidth',1.25)
end
set(gca,'yscale','log')
xlabel('z_i','Interpreter','TeX'); ylabel('\sigma_i(z_i)','Interpreter','TeX'); 
legend('Interpreter','TeX')
ylimtemp = ylim;
ylim([min(10^floor(log10(sigmaSquared(end,1).^.5)), ylimtemp(1)) max(10^ceil(log10(sigmaSquared(1,1).^.5)),ylimtemp(2))])
grid on

subplot(2,1,2); 
bar(sigmaSquared(:,1).^.5)
title('Hankel Singular Values')
xlabel('Order (Number of States)')
ylabel('State Contribution')
set(gca,'yscale','log')
grid on
ylim([10^floor(log10(sigmaSquared(end,1).^.5)) 10^ceil(log10(sigmaSquared(1,1).^.5))])

drawnow

end