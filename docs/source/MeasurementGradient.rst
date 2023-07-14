=========================================================================================================================================
Presentation of the Code that Tests the Model that Uses a Gradient that Determines How Often the Particle has its Wave Function Collapsed
=========================================================================================================================================

As is stated in the title, the code presented below tests the model where the likelihood of the particle located at a particular site having its wave function collapse is dependent on the value of the site's y-index. The higher the value of the y-index (starting at zero), the more likely the particle is to have its wave function collapse. This works by dividing the time evolution of any particular driving step into a series of time steps. After each time step, a y-index is randomly chosen where y-indices that have a higher value are more likely to be chosen than ones that have a lower y-index. Then the eigenvectors and eigenvalues of the relevant density matrix are calculated as well as the probability for the particle to occupy the chosen y-index. A random number generator then decides using the probability calculated if all of the components of the eigenvectors corresponding to the chosen y-index are set to zero or if all of the other components are set to zero. The eigenvectors are then normalized and a new density matrix is reconstructed using the eigenvalues and the new eigenvectors. Finally, the density matrix itself is normalized. The main MATLAB file that runs the sub-files that obtain the data and plots the resulting data is given below:

.. code-block:: matlab

   clear;
   clc;
   % Determine the number of noise realizations that you want to use
   NTot = 100;
   % Determine the number of driving cycles you want to use
   NTime = 1000;
   save('NTime.mat','NTime')
   % Determine the size of the system in the x-direction
   Li = 2;
   % Determine the size of the system in the y-direction
   Lj = 4;
   % Calculate the total number of sites
   numsites = 2+2*(Li-1)+2*Li*(Lj-1);
   % Calculate the number of sites for every value of the y-index
   numj = round(numsites/Lj);
   % Store the probability of the particle occupying each of the j-indices for
   % each driving step and each noise realization for the system unaffected by
   % wave function collapse operations
   jprobs1 = zeros(NTot,NTime,Lj);
   % Store the probability of the particle occupying each of the y-indices for
   % each driving step and each noise realization for the system where
   % wave function collapse operations operations are implemented
   jprobs2 = zeros(NTot,NTime,Lj);
   % Iterate over all of the noise realizations
   for zamp = 1:NTot
       tic
       % Run the file that obtains all of the data
       TwoDimxyQ
       load('jprobsa.mat')
       load('jprobsb.mat')
       % Calculate the probabilities for having the particle occupy each of
       % the y-indices
       for i = 1:NTime
           for j = 1:Lj
               jprobs1(zamp,i,j) = sum(jprobsa(1,((j-1)*numj+1):(j*numj),i));
               jprobs2(zamp,i,j) = sum(jprobsb(1,((j-1)*numj+1):(j*numj),i));
           end
       end
       save('jprobs1.mat','jprobs1')
       save('jprobs2.mat','jprobs2')
       clearvars -except zamp jprobs1 jprobs2 NTot NTime numsites numj Lj
       % Display the time keeping results if you are impatient
       clc
       zamp
       toc
   end
   % Calculate the average probability for the particle to occupy each of the
   % y-indices for each of the driving cycles for the system unaffected by
   % wave function collapse operations
   jprobs1ave = zeros(NTime,Lj);
   % Calculate the average probability for the particle to occupy each of the
   % y-indices for each of the driving cycles for the system where
   % wave function collapse operations are implemented
   jprobs2ave = zeros(NTime,Lj);
   % Calculate the standard deviation associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % unaffected by wave function collapse operations
   jprobs1std = zeros(NTime,Lj);
   % Calculate the standard deviation associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % where wave function collapse operations are implemented
   jprobs2std = zeros(NTime,Lj);
   for i = 1:NTime
       for j = 1:Lj
           jprobs1ave(i,j) = sum(jprobs1(:,i,j))/NTot;
           jprobs2ave(i,j) = sum(jprobs2(:,i,j))/NTot;
           jprobs1std(i,j) = std(jprobs1(:,i,j))/sqrt(NTot-1);
           jprobs2std(i,j) = std(jprobs2(:,i,j))/sqrt(NTot-1);
       end
   end
   % Plot the curves associated with the probability of the particle to occupy
   % each of the y-indices for both the systems that are affected and
   % unaffected by the wave function collapse operations
   for i = 1:Lj
       figure('units','normalized','outerposition',[0 0 1 1]);
       errorbar(1:NTime,jprobs1ave(:,i),jprobs1std(:,i),'Color','b')
       hold on
       errorbar(1:NTime,jprobs2ave(:,i),jprobs2std(:,i),'Color','g')
       hold off
       title(['Probability of Occupying J-index ' num2str(i)],'FontSize',40,'Interpreter','latex')
   end
