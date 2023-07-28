======================================================================================================================
Presentation of the Code that Tests the Model that Uses a Gradient that Determines How Often the Particle is Entangled
======================================================================================================================

As is stated in the title, the code presented below tests the model where the likelihood of the particle located at a particular site becoming entangled with an external qubit is dependent on the value of the site's y-index. The higher the value of the y-index (starting at zero), the more likely the particle is to become entangled with the external qubit. This works by dividing the time evolution of any particular driving step into a series of time steps. After each time step, a y-index is randomly chosen where y-indices that have a higher value are more likely to be chosen than ones that have a lower y-index. Then the algorithm iterates over all possible values of the x-index and :math:`$\alpha$` for the given y-index and adds an external qubit for each iteration. For each iteration, if a particle is present at the site of interest, the external qubit is flipped from spin down to spin up, otherwise it is left alone. Then the external qubit is effectively removed by calculating the reduced density matrix. The main MATLAB file that runs the sub-files that obtain the data and plots the resulting data is given below:

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
   % Calculate the standard error associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % unaffected by entanglement operations
   jprobs1std = zeros(NTime,4);
   % Calculate the standard error associated with the particle occupying
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

The code above runs TwoDimxyQ.m, which is the main file that actually runs the simulation for each noise realization. This code is presented below:

.. code-block:: matlab

    % Define the size of the system
    Li = 2;
    Lj = 4;
    LSquared = 2*Li*Lj;
    % Determine how many qubits are needed to define this system
    nqubits = log2(LSquared);
    % Determine the frequency with which wave function collapse occurs for
    % y-indices 0, 1, 2, and 3
    entprob = [1/10^3 1/10^2 1/10 1];
    % Determine the number of times per driving step that the presence of a
    % particle is measured for a single site
    measint = 100;
    % The following if else statements determines how the time evolution takes
    % place
    if (measint<1)
        timeinterupt = '0';
    else
        timeinterupt = '1';
    end
    % Determine the number of particles that make up the system
    ntimes = 1;
    % Determine the size of the chemical potential
    del = 0.4;
    % Determine the size of the temporal disorder
    tchaos = 0.5;
    % Determine the energy needed for particles to hop between sites
    J = 1.25;
    % NTime sets how many driving cycles the system is evolved for
    load('NTime.mat')
    NVec = 1:NTime;
    N = max(NVec);
    rng('shuffle');
    % The following generates the Hamiltonians for each of the five driving
    % steps
    [H1, H2, H3, H4, H5, V1, V3] = FastTwoDxyHamiltonians(Li,Lj,J,del);
    % Set up the wave function
    W = eye(LSquared);
    wave = W(:,1:ntimes);
    rng('shuffle');
    % Set up the temporal disorder variables for all of the driving steps
    TimeDisorder1 = -tchaos + 2*tchaos*rand(1,N);
    TimeDisorder2 = -tchaos + 2*tchaos*rand(1,N);
    TimeDisorder3 = -tchaos + 2*tchaos*rand(1,N);
    TimeDisorder4 = -tchaos + 2*tchaos*rand(1,N);
    TimeDisorder5 = -tchaos + 2*tchaos*rand(1,N);
    TimeDisorder1 = [-1 TimeDisorder1];
    TimeDisorder2 = [-1 TimeDisorder2];
    TimeDisorder3 = [-1 TimeDisorder3];
    TimeDisorder4 = [-1 TimeDisorder4];
    TimeDisorder5 = [-1 TimeDisorder5];
    wave2 = wave;
    % The following matrix stores information regarding the probability of the
    % particle occupying a site within a particular y-index for the reference
    % system unaffected by entanglement operations
    jprobsa = zeros(1,Lj,N);
    % The following matrix stores information regarding the probability of the
    % particle occupying a site within a particular y-index for the reference
    % system where entanglement operations are implemented
    jprobsb = zeros(1,Lj,N);
    aph = 0;
    % The following matrix stores all of the control operations that flip the
    % external qubit if a qubit is present at a particular site
    measmats = zeros(2^(ntimes*nqubits+1),2^(ntimes*nqubits+1),2*Li*Lj);
    % The following vector stores where the site that becomes entangled with
    % the external qubit changes its y-index value.
    numvec = [];
    for j = (Lj-1):(-1):0
        for i = 0:(Li-1)
            aph = aph + 1;
            % locmat stores the A site of interest for the current iteration of
            % j and i
            locmat = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits));
            % notlocmat stores every site other than the A site of interest for
            % the current iteration of j and i
            notlocmat = eye(2^(ntimes*nqubits),2^(ntimes*nqubits));
            locmat(1+2*i+2*Li*j,1+2*i+2*Li*j) = 1;
            notlocmat(1+2*i+2*Li*j,1+2*i+2*Li*j) = 0;
            % If the particle is present at the A site of interest, flip the
            % external qubit, otherwise leave the external qubit alone.
            measmats(:,:,1+2*i+2*Li*j) = measmats(:,:,1+2*i+2*Li*j) + kron(locmat,[0 1; 1 0]) + kron(notlocmat,[1 0; 0 1]);
            aph = aph + 1;
            % locmat stores the B site of interest for the current iteration of
            % j and i
            locmat = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits));
            % notlocmat stores every site other than the B site of interest for
            % the current iteration of j and i
            notlocmat = eye(2^(ntimes*nqubits),2^(ntimes*nqubits));
            locmat(2+2*i+2*Li*j,2+2*i+2*Li*j) = 1;
            notlocmat(2+2*i+2*Li*j,2+2*i+2*Li*j) = 0;
            % If the particle is present at the B site of interest, flip the
            % external qubit, otherwise leave the external qubit alone.
            measmats(:,:,2+2*i+2*Li*j) = measmats(:,:,2+2*i+2*Li*j) + kron(locmat,[0 1; 1 0]) + kron(notlocmat,[1 0; 0 1]);
        end
        numvec = [numvec aph];
    end
    % Stores how many sites are in the system
    num = aph;
    % Time evolve the system that is unaffected by entanglement operations
    for z = 1:N
        wave2 = expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/5)*expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/5)*expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/5)*expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/5)*expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/5)*wave2;
        % Calculate the probability for the particle to occupy each particular
        % value of the y-index
        for j = 0:(Lj-1)
            probnow = 0;
            for i = 0:(Li-1)
                for k = 1:2
                    probnow = probnow + abs(wave2(k+2*i+2*Li*j))^2;
                end
            end
            jprobsa(1,j+1,z) = probnow;
        end
    end
    % Generate the density matrix for the system where entanglement operations
    % are involved.
    if (ntimes==1)
        density = wave(:,1)*ctranspose(wave(:,1));
    else
        density = kron(wave(:,1)*ctranspose(wave(:,1)),wave(:,2)*ctranspose(wave(:,2)));
        for i = 3:ntimes
            density = kron(density,wave(:,i)*ctranspose(wave(:,i)));
        end
    end
    if (timeinterupt=='1')
        for z = 1:N
            unitnow = expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/(5*measint)));
            end
            for t = 1:measint
                density = unitnow*density*ctranspose(unitnow);
                % Draw a random number to determine probabilities
                draw = rand;
                for t2 = 1:length(entprob)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<entprob(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Iterate over all sites that have the y-index of interest and
                % for all of those sites, add an external qubit. Then, if the
                % particle is present at the site of interest, flip the
                % external qubit from the spin down state to the spin up state.
                % Finally, remove the external qubit by calculating the reduced
                % density matrix.
                for ti = 0:(Li-1)
                    for tk = 1:2
                        density = kron(density,[1 0; 0 0]);
                        density = measmats(:,:,tk+2*ti+2*Li*(cnow-1))*density*ctranspose(measmats(:,:,tk+2*ti+2*Li*(cnow-1)));
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            unitnow = expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/(5*measint)));
            end
            for t = 1:measint
                density = unitnow*density*ctranspose(unitnow);
                % Draw a random number to determine probabilities
                draw = rand;
                for t2 = 1:length(entprob)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<entprob(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Iterate over all sites that have the y-index of interest and
                % for all of those sites, add an external qubit. Then, if the
                % particle is present at the site of interest, flip the
                % external qubit from the spin down state to the spin up state.
                % Finally, remove the external qubit by calculating the reduced
                % density matrix.
                for ti = 0:(Li-1)
                    for tk = 1:2
                        density = kron(density,[1 0; 0 0]);
                        density = measmats(:,:,tk+2*ti+2*Li*(cnow-1))*density*ctranspose(measmats(:,:,tk+2*ti+2*Li*(cnow-1)));
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            unitnow = expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/(5*measint)));
            end
            for t = 1:measint
                density = unitnow*density*ctranspose(unitnow);
                % Draw a random number to determine probabilities
                draw = rand;
                for t2 = 1:length(entprob)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<entprob(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Iterate over all sites that have the y-index of interest and
                % for all of those sites, add an external qubit. Then, if the
                % particle is present at the site of interest, flip the
                % external qubit from the spin down state to the spin up state.
                % Finally, remove the external qubit by calculating the reduced
                % density matrix.
                for ti = 0:(Li-1)
                    for tk = 1:2
                        density = kron(density,[1 0; 0 0]);
                        density = measmats(:,:,tk+2*ti+2*Li*(cnow-1))*density*ctranspose(measmats(:,:,tk+2*ti+2*Li*(cnow-1)));
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            unitnow = expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/(5*measint)));
            end
            for t = 1:measint
                density = unitnow*density*ctranspose(unitnow);
                % Draw a random number to determine probabilities
                draw = rand;
                for t2 = 1:length(entprob)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<entprob(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Iterate over all sites that have the y-index of interest and
                % for all of those sites, add an external qubit. Then, if the
                % particle is present at the site of interest, flip the
                % external qubit from the spin down state to the spin up state.
                % Finally, remove the external qubit by calculating the reduced
                % density matrix.
                for ti = 0:(Li-1)
                    for tk = 1:2
                        density = kron(density,[1 0; 0 0]);
                        density = measmats(:,:,tk+2*ti+2*Li*(cnow-1))*density*ctranspose(measmats(:,:,tk+2*ti+2*Li*(cnow-1)));
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            unitnow = expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/(5*measint)));
            end
            for t = 1:measint
                density = unitnow*density*ctranspose(unitnow);
                % Draw a random number to determine probabilities
                draw = rand;
                for t2 = 1:length(entprob)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<entprob(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Iterate over all sites that have the y-index of interest and
                % for all of those sites, add an external qubit. Then, if the
                % particle is present at the site of interest, flip the
                % external qubit from the spin down state to the spin up state.
                % Finally, remove the external qubit by calculating the reduced
                % density matrix.
                for ti = 0:(Li-1)
                    for tk = 1:2
                        density = kron(density,[1 0; 0 0]);
                        density = measmats(:,:,tk+2*ti+2*Li*(cnow-1))*density*ctranspose(measmats(:,:,tk+2*ti+2*Li*(cnow-1)));
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            % Calculate the probability for the particle to occupy each of the
            % sites
            for j = 0:(Lj-1)
                probnow = 0;
                for i = 0:(Li-1)
                    for k = 1:2
                        probnow = probnow + abs(density(k+2*i+2*Li*j,k+2*i+2*Li*j));
                    end
                end
                jprobsb(1,j+1,z) = probnow;
            end
        end
    else
        % Determine after how many driving steps the entanglement operations
        % are implemented
        measint2 = round(1/measint);
        aph = 0;
        % Iterate over all of the driving cycles
        for z = 1:N
            % Iterate over all of the driving steps
            for z2 = 1:5
                aph = aph + 1;
                % If z2==1, implement the first driving step
                if (z2==1)
                    unitnow = expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % If z2==2, implement the second driving step
                elseif (z2==2)
                    unitnow = expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % If z2==3, implement the third driving step
                elseif (z2==3)
                    unitnow = expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % If z2==4, implement the fourth driving step
                elseif (z2==4)
                    unitnow = expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % If z2==5, implement the fifth driving step
                elseif (z2==5)
                    unitnow = expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                end
                % Draw a random number to determine probabilities
                draw = rand;
                for t2 = 1:length(entprob)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<entprob(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Iterate over all sites that have the y-index of interest and
                % for all of those sites, add an external qubit. Then, if the
                % particle is present at the site of interest, flip the
                % external qubit from the spin down state to the spin up state.
                % Finally, remove the external qubit by calculating the reduced
                % density matrix.
                for ti = 0:(Li-1)
                    for tk = 1:2
                        density = kron(density,[1 0; 0 0]);
                        density = measmats(:,:,tk+2*ti+2*Li*(cnow-1))*density*ctranspose(measmats(:,:,tk+2*ti+2*Li*(cnow-1)));
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
                % If the current iteration is for the fifth driving step,
                % calculate the probability for the particle to occupy each of
                % the y-indices.
                if (z2==5)
                    for j = 0:(Lj-1)
                        probnow = 0;
                        for i = 0:(Li-1)
                            for k = 1:2
                                probnow = probnow + abs(density(k+2*i+2*Li*j,k+2*i+2*Li*j));
                            end
                        end
                        jprobsb(1,j+1,z) = probnow;
                    end
                end
            end
        end
    end
    save('jprobsa.mat','jprobsa')
    save('jprobsb.mat','jprobsb')

This uses the function FastTwoDxyHamiltonians.m, which generates the Hamiltonians that implement the five driving steps. This function is presented as follows:

.. code-block:: matlab

    function [Ham1, Ham2, Ham3, Ham4, Ham5, Vel1, Vel3] = FastTwoDxyHamiltonians(Li,Lj,J,del)
    % This function generates the Hamiltonians that implement the five step
    % Floquet drive as well as the velocity matrices that are used to measure
    % the topological current during the first and third driving steps. The
    % system is defined by Li sites in the x-direction and Lj sites in the
    % y-direction, the hopping strength is given by J, and the strength of the
    % on-site potential implemented during step 5 is given by del.
    %%%
    % Define the total number of sites that defines the system with LSquared
    LSquared = 2*Li*Lj;
    % Initialize all of the Hamiltonians and the velocity matrices as matrices
    % of zeros
    Muy = zeros(LSquared);
    H1 = Muy;
    H2 = Muy;
    H3 = Muy;
    H4 = Muy;
    H5 = Muy;
    V1 = Muy;
    V3 = Muy;
    % Populate all of the Hamiltonians and the velocity matrices in the
    % appropriate locations such that they perform that actions they were
    % intended to.
    for i = 2:2:LSquared
        H1(i,(i-1)) = -J;
        H1((i-1),i) = -J;
        V1((i-1),i) = -1i*J;
        V1(i,(i-1)) = 1i*J;
    end
    clear i
    for i = 0:(Li-1)
        for j = 0:(Lj-2)
            H2((2+2*i+2*Li*(j+1)),(1+2*rem((i+1),Li)+2*Li*j)) = -J;
            H2((1+2*rem((i+1),Li)+2*Li*j),(2+2*i+2*Li*(j+1))) = -J;
            H4((2+2*i+2*Li*j),(1+2*i+2*Li*(j+1))) = -J;
            H4((1+2*i+2*Li*(j+1)),(2+2*i+2*Li*j)) = -J;
        end
        clear j
        for j = 0:(Lj-1)
            H3((1+2*rem((i+1),Li)+2*Li*j),(2+2*i+2*Li*j)) = -J;
            H3((2+2*i+2*Li*j),(1+2*rem((i+1),Li)+2*Li*j)) = -J;
            V3((1+2*rem((i+1),Li)+2*Li*j),(2+2*i+2*Li*j)) = -1i*J;
            V3((2+2*i+2*Li*j),(1+2*rem((i+1),Li)+2*Li*j)) = 1i*J;
        end
    end
    for k = 1:LSquared
        H5(k,k) = ((-1)^(k-1))*del;
    end
    % Give the results as output.
    Ham1 = H1;
    Ham2 = H2;
    Ham3 = H3;
    Ham4 = H4;
    Ham5 = H5;
    Vel1 = V1;
    Vel3 = V3;
    end

An additional helper function named ReducedDensity.m is used to calculate the reduced density matrix and thereby, effectively remove the additional qubit.

.. code-block:: matlab

    function [rdensity] = ReducedDensity(densityi,size,targets)
    % This function takes the density matrix densityi composed of size qubits
    % and calculates the reduced density matrix for the qubits given by targets
    % and returns this reduced density matrix as rdensity
    %%%
    % Determine the number of qubits that compose targets
    nq = length(targets);
    % Determine the number of qubits in densityi that are not going to compose
    % the outputted reduced density matrix
    nq2 = size - nq;
    % Initialize the matrix that will store the reduced density matrix
    redden = zeros(2^nq);
    % Iterate over all possible configurations of the qubits that will not
    % compose the reduced density matrix
    for i = 0:(2^nq2-1)
        % Express the number for the current iteration as a bitstring of length
        % nq2
        const = dec2bin(i);
        const2 = nq2 - length(const);
        for j = 1:const2
            const = ['0' const];
        end
        % count is used to determine how far across the bitstring we have gone
        % when using the information in the bitstring to generate the matrix
        % opmat that will be used to create the reduced density matrix.
        count = 0;
        % If 1 is an entry of targets, then make the first matrix that composes
        % the set of Kronecker products that generates opmat be the 2 by 2
        % identity matrix
        if sum(1==targets)
            opmat = eye(2);
        else
        % Otherwise make the first matrix that composes this set of Kronecker
        % products be the appropriate single qubit spin vector
            count = count+1;
            if (const(count)=='1')
                opmat = [0; 1];
            else
                opmat = [1; 0];
            end
        end
        % Iterate through all of the rest of the qubits (both the target qubits
        % for the reduced density matrix as well as all of the other qubits)
        % and determine whether the next matrix in the set of Kronecker
        % products should be an identity matrix or the spin up or down state
        % vector. If the qubit of interest is a target qubit for the reduced
        % density matrix then use the identity matrix otherwise use the
        % appropriate state vector.
        for j = 2:size
            if sum(j==targets)
                opmat = kron(opmat,eye(2));
            else
                count = count + 1;
                if (const(count)=='1')
                    opmat = kron(opmat,[0; 1]);
                else
                    opmat = kron(opmat,[1; 0]);
                end
            end
        end
        % Use opmat to perform operations on densityi in order to obtain the
        % appropriate information about the reduced density matrix and add this
        % information to redden.
        redden = redden + ctranspose(opmat)*densityi*opmat;
    end
    % Normalize redden
    redden = redden/trace(abs(redden));
    % Return the reduced density matrix as rdensity
    rdensity = redden;
    end
