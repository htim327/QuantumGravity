==============================================================================================================================================================
Presentation of the Code that Tests the Model that Uses a Gradient that Determines How Often the Particle is both Entangled and has its Wave Function Collapse
==============================================================================================================================================================

The following code presents a method where the time evolution of each driving step is divided into a series of time steps. After each time step, the algorithm undergoes a series of iterations where during each iteration, the particle is entangled with an external qubit if the particle is present at that site and then the particle has its wave function collapse to a single site if the particle is located at a site site of interest. For both of these process, a random value of the y-index is obtained (separately for each process) where a y-index is more likely to be chosen if it has a higher value with the minimum value being zero. In addition, both of these processes separately calculate a random value for the x-index and :math:`$\alpha$` using a uniform distribution. The number of iterations was arbitrarily chosen to be equal to the number of sites of the system. For the entanglement process, first an external qubit is added. Then if the particle is present at the site of interest, the external qubit is flipped from spin down to spin up. Finally, the particle is effectively removed by calculating the reduced density matrix. For the wave function collapse process, first the eigenvalues and eigenvectors of the relevant density matrix are calculated and the probability for the particle to occupy the site of interest is calculated. Then a random number generator is used along with the probability calculated to determine whether to set the entries corresponding to the site of interest to zero or to set all of the other entries to zero. Then a new density matrix is calculated using the eigenvalues and the new eigenvectors and the new density matrix is normalized. The main MATLAB file that runs the sub-files that obtain the data and plots the resulting data is given below:

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
   % entanglement and wave function collapse
   jprobs1 = zeros(NTot,NTime,Lj);
   % Store the probability of the particle occupying each of the y-indices for
   % each driving step and each noise realization for the system where
   % entanglement and wave function collapse operations are implemented
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
       % Display the timekeeping results if you are impatient
       clc
       zamp
       toc
   end
   % Calculate the average probability for the particle to occupy each of the
   % y-indices for each of the driving cycles for the system unaffected by
   % entanglement and wave function collapse
   jprobs1ave = zeros(NTime,Lj);
   % Calculate the average probability for the particle to occupy each of the
   % y-indices for each of the driving cycles for the system where
   % entanglement and wave function collapse operations are implemented
   jprobs2ave = zeros(NTime,Lj);
   % Calculate the standard deviation associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % unaffected by entanglement and wave function collapse
   jprobs1std = zeros(NTime,Lj);
   % Calculate the standard deviation associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % where entanglement and wave function collapse operations are implemented
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
   % unaffected by the entanglement operations
   for i = 1:Lj
       figure('units','normalized','outerposition',[0 0 1 1]);
       errorbar(1:NTime,jprobs1ave(:,i),jprobs1std(:,i),'Color','b')
       hold on
       errorbar(1:NTime,jprobs2ave(:,i),jprobs2std(:,i),'Color','g')
       hold off
       title(['Probability of Occupying J-index ' num2str(i)],'FontSize',40,'Interpreter','latex')
   end
