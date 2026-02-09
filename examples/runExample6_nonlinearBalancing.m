function runExample6_nonlinearBalancing()
%runExample6_nonlinearBalancing Runs the finite element beam example
%to compare linear vs nonlinear balancing approaches.
%
%   Usage:  runExample6_nonlinearBalancing()
%
%   Description: The cantilever Euler-Bernoulli beam with von Karman geometric
%   nonlinearity was first used as an example in the context of nonlinear
%   balancing [1], where we were just using it as a scalable example with
%   polynomial nonlinearities to illustrate computing energy functions. Here, we
%   come full circle and consider the problem for computing balanced
%   realizations and ultimately reduced-order models.
%
%   This problem is interesting because the nonlinear model is not just a smooth
%   slightly perturbed version of the linear model, as is the case for example
%   with the pendulum. In this case, there are qualitative features of the
%   nonlinear model which are not captured by the linearized dynamics, not just
%   due to a nonlinear output equation, but rather due to nonlinear couplings in
%   the dynamics. Specifically, the linearized beam theory cannot exhibit
%   coupling between transverse and axial deformations, whereas the nonlinear
%   theory does. The first plot produced illustrates this.
%
%   The next set of results is produced with a discretization consisting of a
%   single element, for which the governing equations are of order 6. This is
%   sufficiently small that we can actually inspect the governing equations
%   directly and observe some details of what the balancing transformation is
%   doing, and what is omitted in a linear balancing transformation with regard
%   to the nonlinear coupling in the dynamics. We show that the balancing
%   transformation performance is dictated by two key issues: 1) bijectivity of
%   the balancing transformation, and 2) the ability to encode the coupling
%   between states in order to eliminate redundancy. Linear balancing
%   transformations achieve bijectivity more straightforwardly than nonlinear
%   transformations, but this problem illustrates the inability to capture some
%   kinds of nonlinear couplings in the dynamics. As a result, once truncation
%   is considered as a means of producing a reduced-order model, the linear
%   balancing procedure performs worse locally. This is visible directly in the
%   governing equations.
%
%   We also illustrate this characteristic on simulation results for the
%   reduced-order models obtained by linear balancing and nonlinear balancing.
%   Our hypothesis is that the nonlinear balancing transformations more
%   accurately capture the coupling and redundancy in the dynamics, so locally
%   higher-order transformations should perform better. However, the nonlinear
%   balancing transformations are only convergent locally, so as we consider
%   larger initial conditions, we expect to see the nonlinear solutions diverge
%   eventually, which we show through a series of three different initial
%   conditions of increasing size.
%
%   We compute the energy functions, the input-normal/output-diagonal
%   transformation, and then the true balancing transformation, given by the
%   composition x = ̅Φ(z̄(z̄) = Φ(𝝋(z̄)). We visualize this mapping
%   from the z̄ coordinates to the x coordinates by forming a grid in the
%   z̄ coordinates and mapping that grid to the x coordinates. This is
%   compared for linear and nonlinear transformations, and the polynomial
%   approach is compared with the more direct Newton iteration approach to
%   illustrate the similarity in terms of accuracy but major speedup with
%   the polynomial approach.
%
%   References: [1]
%
%   Part of the NLbalancing repository.
%%

% %% N=1 elements (n=6)
% % First, we show that the linearized dynamics are qualitatively different
% % than the nonlinear dynamics
% runExample6_linearComparison()
% 
% % Next, we show the differences in the models obtained via linear and nonlinear
% % balancing. These are all nonlinear models; the difference here is specifically
% % in the degree of balancing transformation used.
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,6,1e3,true,12,plot=true);
% runExample6_balancedReduction(x0,3,1,6,1e3,true,12,plot=true);
% runExample6_balancedReduction(x0,4,1,6,1e3,true,12,plot=true);
% 
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,4,1e3,false,13,plot=true);
% runExample6_balancedReduction(x0,3,1,4,1e3,false,13,plot=true);
% runExample6_balancedReduction(x0,4,1,4,1e3,false,13,plot=true);
% 
% % Now we show that, as the initial conditions depart from the origin, eventually
% % the nonlinear approximations are no longer convergent.
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,6,1e4,false,14,plot=true);
% runExample6_balancedReduction(x0,3,1,6,1e4,false,14,plot=true);
% runExample6_balancedReduction(x0,4,1,6,1e4,false,14,plot=true);
% 
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,4,1e4,false,15,plot=true);
% runExample6_balancedReduction(x0,3,1,4,1e4,false,15,plot=true);
% runExample6_balancedReduction(x0,4,1,4,1e4,false,15,plot=true);
% 
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,6,2e4,false,16,plot=true);
% runExample6_balancedReduction(x0,3,1,6,2e4,false,16,plot=true);
% runExample6_balancedReduction(x0,4,1,6,2e4,false,16,plot=true);
% 
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,4,2e4,false,17,plot=true);
% runExample6_balancedReduction(x0,3,1,4,2e4,false,17,plot=true);
% runExample6_balancedReduction(x0,4,1,4,2e4,false,17,plot=true);
% 
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,6,3e4,false,18,plot=true);
% runExample6_balancedReduction(x0,3,1,6,3e4,false,18,plot=true);
% runExample6_balancedReduction(x0,4,1,6,3e4,false,18,plot=true);
% 
% [~,~,~,x0] = runExample6_balancedReduction([],2,1,4,3e4,false,19,plot=true);
% runExample6_balancedReduction(x0,3,1,4,3e4,false,19,plot=true);
% runExample6_balancedReduction(x0,4,1,4,3e4,false,19,plot=true);

%% N=N elements (scalability results)
% Next we compare the time required to simulate the FOM (T2), construct the ROM
% (T1), and simulate the ROM (T3).
U0 = 2e4;
numEls = [1 2 4 8 16 32 64 128 180];
x0s = runExample6_getStaticDeflectionIC(numEls, U0);

Linear balancing transformation
runExample6_timeTrials(U0, numEls, 2, 9, x0s);

% Quadratic balancing transformation
runExample6_timeTrials(U0, numEls, 3, 9, x0s);

% Cubic balancing transformation
runExample6_timeTrials(U0, numEls, 4, 8, x0s);

end

function runExample6_timeTrials(U0, numEls, d, numTrials, x0s)
Ex6timings = zeros(numTrials,3);
for i = 1:min([3,numTrials])
    T1Temp = 0; T2Temp = 0; T3Temp = 0;
    for j=1:3 % average over 3 runs
        [T1, T2, T3] = runExample6_balancedReduction(x0s{i},d,numEls(i),6,U0,false,1,plot=false);
        T1Temp = T1Temp+T1; T2Temp = T2Temp+T2; T3Temp = T3Temp+T3;
    end
    Ex6timings(i,:) = [T1Temp, T2Temp, T3Temp]./3;
end
for i = 4:numTrials
    try
        zeros((6*numEls(i))^d,1);
        [T1, T2, T3] = runExample6_balancedReduction(x0s{i},d,numEls(i),6,U0,false,1,plot=false);
        Ex6timings(i,:) = [T1, T2, T3];
    catch
        warning('RAM capacity will be exceeded, skipping this case')
    end
end
fprintf('Writing data to plots/example6_balancingScaling_d%d.dat \n',d)
fileID = fopen(sprintf('plots/example6_balancingScaling_d%d.dat',d), 'w');
fprintf(fileID, '# Table I Data\n# finite element beam model, scalability results; d=%d \nnumElements &   n  &     n^%d      &  Balancing CPU-sec   &  FOM Sim. CPU-sec    &  ROM Sim. CPU-sec \n', d, d+1);
for i=1:numTrials
    fprintf(fileID, ' %5d      &%4d  &  %10.4e  &     %12.6e     &     %12.6e     &    %12.6e \n', numEls(i), 6*numEls(i), (6*numEls(i))^(d+1), Ex6timings(i,1), Ex6timings(i,2), Ex6timings(i,3));
end
fclose(fileID);
end