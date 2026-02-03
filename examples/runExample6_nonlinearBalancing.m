function runExample6_nonlinearBalancing()

% N=1 elements (n=6):
% - first, show that the linearization of these dynamics is qualitatively
% different than the nonlinear dynamics at this initial condition
runExample6_linearComparison()

% - print output equation to command window containing nonlinearities that
% become truncated
% - small initial condition, quadratic and cubic transformations similar
[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,6,1e3,true,12,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,6,1e3,true,12,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,6,1e3,true,12,plot=true);

[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,4,1e3,false,13,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,4,1e3,false,13,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,4,1e3,false,13,plot=true);

% - larger initial condition, cubic transformations slightly better than
% quadratic as they both break down, and the full transformation itself 
% begins to no-longer be bijective at these points in the state space
[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,6,1e4,false,14,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,6,1e4,false,14,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,6,1e4,false,14,plot=true);

[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,4,1e4,false,15,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,4,1e4,false,15,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,4,1e4,false,15,plot=true);

[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,6,2e4,false,16,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,6,2e4,false,16,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,6,2e4,false,16,plot=true);

[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,4,2e4,false,17,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,4,2e4,false,17,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,4,2e4,false,17,plot=true);

[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,6,3e4,false,18,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,6,3e4,false,18,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,6,3e4,false,18,plot=true);

[~,~,~,~,x0] = runExample6_balancedReduction_staticDeflectionIC([],3,2,1,4,3e4,false,19,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,3,1,4,3e4,false,19,plot=true);
runExample6_balancedReduction_staticDeflectionIC(x0,3,4,1,4,3e4,false,19,plot=true);


% N=3 elements (n=18):
% - output equation containing nonlinearities that become truncated
% - small initial condition, quadratic and cubic transformations similar
runExample6_balancedReduction_staticDeflectionIC([],3,2,3,18,1e3,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,3,3,18,1e3,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,4,3,18,1e3,false,plot=true);

runExample6_balancedReduction_staticDeflectionIC([],3,2,3,14,1e3,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,3,3,14,1e3,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,4,3,14,1e3,false,plot=true);

% - larger initial condition, cubic transformations slightly better than
% quadratic as they both break down, and the full transformation itself 
% begins to no-longer be bijective at these points in the state space
runExample6_balancedReduction_staticDeflectionIC([],3,2,3,14,1e4,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,3,3,14,1e4,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,4,3,14,1e4,false,plot=true);

% instability prone (these fail)
% *** The reason I've been having such a hard time here is its the initial
% condition that ends up causing the instability; this can be seen in the
% N=1 element sims, so perhaps that's best to explain there somehow
runExample6_balancedReduction_staticDeflectionIC([],3,2,3,14,9e4,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,3,3,14,9e4,false,plot=true);
runExample6_balancedReduction_staticDeflectionIC([],3,4,3,14,9e4,false,plot=true);

% Instead of making the error plot using the L2 output error, can I just
% look at the error of the initial condition (step output of FOM) through
% the balancing transformation? I.e. the error of passing it to and back





% N=1 elements (n=6) to show qualitatively the output equation containing nonlinearities that become truncated
runExample6_balancedReduction(2,1,4,1e-3,true,plot=true);
runExample6_balancedReduction(3,1,4,1e-3,true,plot=true);
runExample6_balancedReduction(4,1,4,1e-3,true,plot=true);

% N=3 (r=18) elements, different r; garbage in, garbage out: if the FOM
% transformed is no good, then the ROM is even worse! This is dependent on
% the distance to the origin.

% For this state region, the transformation is valid -> good ROM result
U0 = 1e-3;
runExample6_balancedReduction(2,3,18,U0,false,plot=true);
runExample6_balancedReduction(3,3,18,U0,false,plot=true);
runExample6_balancedReduction(4,3,18,U0,false,plot=true);

runExample6_balancedReduction(2,3,16,U0,false,plot=true);
runExample6_balancedReduction(3,3,16,U0,false,plot=true);
runExample6_balancedReduction(4,3,16,U0,false,plot=true);

runExample6_balancedReduction(2,3,6,U0,false,plot=true);
runExample6_balancedReduction(3,3,6,U0,false,plot=true);
runExample6_balancedReduction(4,3,6,U0,false,plot=true);
% For this state region, the transformation is invalid -> bad ROM result
U0 = 5e-1;
runExample6_balancedReduction(2,3,18,U0,false,plot=true);
runExample6_balancedReduction(3,3,18,U0,false,plot=true);
runExample6_balancedReduction(4,3,18,U0,false,plot=true);

runExample6_balancedReduction(2,3,16,U0,false,plot=true);
runExample6_balancedReduction(3,3,16,U0,false,plot=true);
runExample6_balancedReduction(4,3,16,U0,false,plot=true);

runExample6_balancedReduction(2,3,6,U0,false,plot=true);
runExample6_balancedReduction(3,3,6,U0,false,plot=true);
runExample6_balancedReduction(4,3,6,U0,false,plot=true);


% Illustration that error is lower locally but higher far away 
degrees = [2 3 4]; U0s = linspace(1e-3,5e-2,20); r = 18;
errors = NaN(length(degrees),length(U0s));
for i = 1:length(degrees)
    for j = 1:length(U0s)
        try
            [~,~,~,errors(i,j)] = runExample6_balancedReduction(degrees(i),3,r,U0s(j),false,"plot",false);
        catch
        end
    end
end
figure; hold on;
for i = 1:length(degrees)
    plot(U0s,errors(i,:),'DisplayName',sprintf('Degree %i transformation',degrees(i)-1))
end
legend


degrees = [2 3 4]; U0s = linspace(1e-3,1e-2,20); r = 16;
errors = NaN(length(degrees),length(U0s));
for i = 1:length(degrees)
    for j = 1:length(U0s)
        try
            [~,~,~,errors(i,j)] = runExample6_balancedReduction(degrees(i),3,r,U0s(j),false,"plot",false);
        catch
        end
    end
end
figure; hold on;
for i = 1:length(degrees)
    plot(U0s,errors(i,:),'DisplayName',sprintf('Degree %i transformation',degrees(i)-1))
end
legend


degrees = [2 3 4]; U0s = linspace(1e-3,5e-2,20); r = 6;
errors = NaN(length(degrees),length(U0s));
for i = 1:length(degrees)
    for j = 1:length(U0s)
        try
            [~,~,~,errors(i,j)] = runExample6_balancedReduction(degrees(i),3,r,U0s(j),false,"plot",false);
        catch
        end
    end
end
figure; hold on;
for i = 1:length(degrees)
    plot(U0s,errors(i,:),'DisplayName',sprintf('Degree %i transformation',degrees(i)-1))
end
legend



runExample6_balancedReduction(2,3,16,1e-1,false,plot=true);
runExample6_balancedReduction(3,3,16,1e-1,false,plot=true);
runExample6_balancedReduction(4,3,16,1e-1,false,plot=true);


runExample6_balancedReduction(2,3,6,4e-1,false,plot=true);
runExample6_balancedReduction(3,3,6,4e-1,false,plot=true);
runExample6_balancedReduction(4,3,6,4e-1,false,plot=true);

% Illustration that error is lower locally but higher far away 
degrees = [2 3 4]; U0s = linspace(1e-3,4e-1,20);
errors = NaN(length(degrees),length(U0s));
for i = 1:length(degrees)
    for j = 1:length(U0s)
        try
            [~,~,~,errors(i,j)] = runExample6_balancedReduction(degrees(i),3,16,U0s(j),false,"plot",false);
        catch
        end
    end
end
figure; hold on;
for i = 1:length(degrees)
    plot(U0s,errors(i,:),'DisplayName',sprintf('Degree %i transformation',degrees(i)-1))
end
legend

runExample6_balancedReduction(2,4,12,1e0);
runExample6_balancedReduction(3,4,12,1e0);

% Illustration that error is lower locally but higher far away 
degrees = [2 3 4]; U0s = linspace(1e-3,1.5e-1,50);
errors = NaN(length(degrees),length(U0s));
for i = 1:length(degrees)
    for j = 1:length(U0s)
        try
            [~,~,~,errors(i,j)] = runExample6_balancedReduction(degrees(i),1,4,U0s(j),false,"plot",false);
        catch
        end
    end
end
figure; hold on;
for i = 1:length(degrees)
    plot(U0s,errors(i,:),'DisplayName',sprintf('Degree %i transformation',degrees(i)-1))
end
legend
runExample6_balancedReduction(2,1,4,U0s(1),false,"plot",true);
runExample6_balancedReduction(3,1,4,U0s(1),false,"plot",true);
runExample6_balancedReduction(4,1,4,U0s(1),false,"plot",true);
runExample6_balancedReduction(2,1,4,U0s(end),false,"plot",true);
runExample6_balancedReduction(3,1,4,U0s(end),false,"plot",true);
runExample6_balancedReduction(4,1,4,U0s(end),false,"plot",true);


% Timing: T1 is the balancing computation, T2 is the FOM simulation time,
% T3 is the ROM simulation time
numEls = [1 2 4 8 16];
Ex6timings = zeros(length(numEls),3);
for i = 1:4
    T1Temp = 0; T2Temp = 0; T3Temp = 0;
    for j=1:3 % average over 3 runs
        [T1, T2, T3] = runExample6_balancedReduction(4, numEls(i), 8);
        T1Temp = T1Temp+T1; T2Temp = T2Temp+T2; T3Temp = T3Temp+T3;
    end
    Ex6timings(i,:) = [T1Temp, T2Temp, T3Temp]./3;
end
for i = 5:length(numEls)
    [T1, T2, T3] = runExample6_balancedReduction(4, numEls(i), 8);
    Ex6timings(i,:) = [T1, T2, T3];
end

fprintf('Writing data to plots/example6_balancingScaling_d4.dat \n')
fileID = fopen('plots/example6_balancingScaling_d4.dat', 'w');
fprintf(fileID, '# Table I Data\n# finite element beam model, scalability results; d=%d \nnumElements &   n  &     n^%d      &  Balancing CPU-sec   &  FOM Sim. CPU-sec    &  ROM Sim. CPU-sec \n', 4, 5);
for i=1:length(numEls)
    fprintf(fileID, ' %5d      &%4d  &  %10.4e  &     %12.6e     &     %12.6e     &    %12.6e \n', numEls(i), 6*numEls(i), (6*numEls(i))^5, Ex6timings(i,1), Ex6timings(i,2), Ex6timings(i,3));
end
fclose(fileID);


numEls = [1 2 4 8 16];
Ex6timings = zeros(length(numEls),3);
for i = 1:4
    T1Temp = 0; T2Temp = 0; T3Temp = 0;
    for j=1:3 % average over 3 runs
        [T1, T2, T3] = runExample6_balancedReduction(3, numEls(i), 8);
        T1Temp = T1Temp+T1; T2Temp = T2Temp+T2; T3Temp = T3Temp+T3;
    end
    Ex6timings(i,:) = [T1Temp, T2Temp, T3Temp]./3;
end
for i = 5:length(numEls)
    [T1, T2, T3] = runExample6_balancedReduction(3, numEls(i), 8);
    Ex6timings(i,:) = [T1, T2, T3];
end

fprintf('Writing data to plots/example6_balancingScaling_d3.dat \n')
fileID = fopen('plots/example6_balancingScaling_d3.dat', 'w');
fprintf(fileID, '# Table I Data\n# finite element beam model, scalability results; d=%d \nnumElements &   n  &     n^%d      &  Balancing CPU-sec   &  FOM Sim. CPU-sec    &  ROM Sim. CPU-sec \n', 3, 4);
for i=1:length(numEls)
    fprintf(fileID, ' %5d      &%4d  &  %10.4e  &     %12.6e     &     %12.6e     &    %12.6e \n', numEls(i), 6*numEls(i), (6*numEls(i))^4, Ex6timings(i,1), Ex6timings(i,2), Ex6timings(i,3));
end
fclose(fileID);

end