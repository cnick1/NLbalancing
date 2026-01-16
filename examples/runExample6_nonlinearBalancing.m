function runExample6_nonlinearBalancing()

% N=1 elements (n=6) to show qualitatively the output equation containing nonlinearities that become truncated
runExample6_balancedReduction(2,1,4,1.5e-2,true);
runExample6_balancedReduction(3,1,4,1.5e-2,true);
runExample6_balancedReduction(4,1,4,1.5e-2,true);

runExample6_balancedReduction(2,4,12,1e0);
runExample6_balancedReduction(3,4,12,1e0);

% Illustration that error is lower locally but higher far away 
degrees = [2 3 4]; U0s = 1:2:100;
errors = zeros(length(degrees),length(U0s));
for i = 1:length(degrees)
    for j = 1:length(U0s)
        [~,~,~,errors(i,j)] = runExample6_balancedReduction(degrees(i),4,12,U0s(j)*1e-2,false);
    end
end
figure; hold on;
for i = 1:length(degrees)
    plot(U0s,errors(i,:),'DisplayName',sprintf('Degree %i transformation',degrees(i)-1))
end
legend
runExample6_balancedReduction(2,4,12,U0s(1)*1e-2,false,"plot",true);
runExample6_balancedReduction(3,4,12,U0s(1)*1e-2,false,"plot",true);
runExample6_balancedReduction(4,4,12,U0s(1)*1e-2,false,"plot",true);
runExample6_balancedReduction(2,4,12,U0s(end)*1e-2,false,"plot",true);
runExample6_balancedReduction(3,4,12,U0s(end)*1e-2,false,"plot",true);
runExample6_balancedReduction(4,4,12,U0s(end)*1e-2,false,"plot",true);


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

fprintf('Writing data to plots/example6_balancingScaling.dat \n')
fileID = fopen('plots/example6_balancingScaling.dat', 'w');
fprintf(fileID, '# Table I Data\n# finite element beam model, scalability results; d=%d \nnumElements &   n  &     n^%d      &  Balancing CPU-sec   &  FOM Sim. CPU-sec    &  ROM Sim. CPU-sec \n', 4, 5);
for i=1:length(numEls)
    fprintf(fileID, ' %5d      &%4d  &  %10.4e  &     %12.6e     &     %12.6e     &    %12.6e \n', numEls(i), 6*numEls(i), (6*numEls(i))^5, Ex6timings(i,1), Ex6timings(i,2), Ex6timings(i,3));
end
fclose(fileID);

end