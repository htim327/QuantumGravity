======================================================================================================================
Presentation of the Code that Tests the Model that Uses a Gradient that Determines How Often the Particle is Entangled
======================================================================================================================

As is stated in the title, the code presented below tests the model where the likelihood of the particle located at a particular site becoming entangled with an external qubit is dependent on the value of the site's y-index. The higher the value of the y-index (starting at zero), the more likely the particle is to become entangled with the external qubit. This works by dividing the time evolution of any particular driving step into a series of time steps. After each time step, a y-index is randomly chosen where y-indices that have a higher value are more likely to be chosen than ones that have a lower y-index. Then the algorithm iterates over all possible values of the x-index and the value for :math:`$\alpha$` for the given y-index and adds an external qubit for each iteration. For each iteration, if a particle is present at the site of interest, the external qubit is flipped from spin down to spin up, otherwise it is left alone. Then the external qubit is left alone. The main MATLAB file that runs the sub-files that obtain the data and plots the resulting data is given below:

.. code-block:: matlab

   clear;
   clc;
   % Determine the number of noise realizations that you want to use
   NTot = 100;
   % Determine the number of driving cycles you want to use
   NTime = 1000;
   save('NTime.mat','NTime')
   % Store the probability of the particle occupying each of the j-indices for
   % each driving step and each noise realization for the system unaffected by
   % entanglement operations
   jprobs1 = zeros(NTot,NTime,4);
   % Store the probability of the particle occupying each of the y-indices for
   % each driving step and each noise realization for the system where
   % entanglement operations are implemented
   jprobs2 = zeros(NTot,NTime,4);
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
           for j = 1:4
               jprobs1(zamp,i,j) = jprobsa(1,j,i);
               jprobs2(zamp,i,j) = jprobsb(1,j,i);
           end
       end
       save('jprobs1.mat','jprobs1')
       save('jprobs2.mat','jprobs2')
       clearvars -except zamp jprobs1 jprobs2 NTot NTime
       % Display the time keeping results if you are impatient
       clc
       zamp
       toc
   end
   % Calculate the average probability for the particle to occupy each of the
   % y-indices for each of the driving cycles for the system unaffected by
   % entanglement operations
   jprobs1ave = zeros(NTime,4);
   % Calculate the average probability for the particle to occupy each of the
   % y-indices for each of the driving cycles for the system where
   % entanglement operations are implemented
   jprobs2ave = zeros(NTime,4);
   % Calculate the standard deviation associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % unaffected by entanglement operations
   jprobs1std = zeros(NTime,4);
   % Calculate the standard deviation associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % where entanglement are implemented
   jprobs2std = zeros(NTime,4);
   for i = 1:NTime
       for j = 1:4
           jprobs1ave(i,j) = sum(jprobs1(:,i,j))/NTot;
           jprobs2ave(i,j) = sum(jprobs2(:,i,j))/NTot;
           jprobs1std(i,j) = std(jprobs1(:,i,j))/sqrt(NTot-1);
           jprobs2std(i,j) = std(jprobs2(:,i,j))/sqrt(NTot-1);
       end
   end
   % Plot the curves associated with the probability of the particle to occupy
   % each of the y-indices for both the systems that are affected and
   % unaffected by the entanglement operations
   for i = 1:4
       figure('units','normalized','outerposition',[0 0 1 1]);
       errorbar(1:NTime,jprobs1ave(:,i),jprobs1std(:,i),'Color','b')
       hold on
       errorbar(1:NTime,jprobs2ave(:,i),jprobs2std(:,i),'Color','g')
       hold off
       title(['Probability of Occupying J-index ' num2str(i)],'FontSize',40,'Interpreter','latex')
   end
