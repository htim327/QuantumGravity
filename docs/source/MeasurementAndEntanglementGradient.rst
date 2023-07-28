==============================================================================================================================================================
Presentation of the Code that Tests the Model that Uses a Gradient that Determines How Often the Particle is both Entangled and has its Wave Function Collapse
==============================================================================================================================================================

The following code presents a method where the time evolution of each driving step is divided into a series of time steps. After each time step, the algorithm undergoes a series of iterations where during each iteration, the particle has its wave function collapse to a single site if the particle is located at a site of interest that is randomly chosen and then the particle is entangled with an external qubit if the particle is present at another randomly chosen site of interest. For both of these process, a random value of the y-index is obtained (separately for each process) where a y-index is more likely to be chosen if it has a higher value with the minimum value being zero. In addition, both of these processes separately calculate a random value for the x-index and :math:`$\alpha$` using a uniform distribution. The number of iterations was arbitrarily chosen to be equal to the number of sites of the system. For the entanglement process, first an external qubit is added. Then if the particle is present at the site of interest, the external qubit is flipped from spin down to spin up. Finally, the particle is effectively removed by calculating the reduced density matrix. For the wave function collapse process, first the eigenvalues and eigenvectors of the relevant density matrix are calculated and the probability for the particle to occupy the site of interest is calculated. Then a random number generator is used along with the probability calculated to determine whether to set the entry corresponding to the site of interest to zero or to set all of the other entries to zero. Then a new density matrix is calculated using the eigenvalues and the new eigenvectors and the new density matrix is also normalized. The main MATLAB file that runs the sub-files that obtain the data and plots the resulting data is given below:

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
   % Calculate the standard error associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % unaffected by entanglement and wave function collapse
   jprobs1std = zeros(NTime,Lj);
   % Calculate the standard error associated with the particle occupying
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
   % unaffected by the entanglement and wave function collapse operations
   for i = 1:Lj
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
    probvec = [1/10^3 1/10^2 1/10 1];
    % Determine the frequency with which entanglement with an external particle
    % occurs for y-indices 0, 1, 2, and 3
    entprob = [1/10^3 1/10^2 1/10 1];
    % Determine the number of times per driving step that a single site is
    % entangled with an external particle and the presence of a particle is
    % measured for a single site
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
    % steps. The velocity matrices V1 and V3 are not useful for this calculation.
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
    % particle occupying each sites for the reference system unaffected by
    % entanglement or wave function collapse
    jprobsa = zeros(1,2+2*(Li-1)+2*Li*(Lj-1),N);
    % The following matrix stores information regarding the probability of the
    % particle occupying each sites for the system where entanglement and wave
    % function collapse are involved
    jprobsb = zeros(1,2+2*(Li-1)+2*Li*(Lj-1),N);
    aph = 0;
    % The following matrix stores the projection operators that are used to
    % calculate the probability of the particle occupying each of the sites
    sitexpectations = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits),2+2*(Li-1)+2*Li*(Lj-1));
    for j = 0:(Lj-1)
        for i = 0:(Li-1)
            for k = 1:2
                aph = aph + 1;
                sitexpectations(k+2*i+2*Li*j,k+2*i+2*Li*j,aph) = 1;
            end
        end
    end
    % Stores how many sites are in the system
    num = aph;
    aph = 0;
    % The following matrix stores all of the control operations that flip the
    % external qubit if a qubit is present at a particular site
    measmats = zeros(2^(ntimes*nqubits+1),2^(ntimes*nqubits+1),2*Li*Lj);
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
    end
    % Time evolve the system that is unaffected by wave function collapse and
    % entanglement
    for z = 1:N
        wave2 = expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/5)*expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/5)*expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/5)*expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/5)*expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/5)*wave2;
        % Calculate the probability for the particle to occupy each of the
        % sites
        for j = 1:num
            jprobsa(1,j,z) = ctranspose(wave2)*sitexpectations(:,:,j)*wave2;
        end
    end
    % Generate the density matrix for the system where both entanglement and
    % wave function collapse are involved.
    if (ntimes==1)
        density = wave(:,1)*ctranspose(wave(:,1));
    else
        density = kron(wave(:,1)*ctranspose(wave(:,1)),wave(:,2)*ctranspose(wave(:,2)));
        for i = 3:ntimes
            density = kron(density,wave(:,i)*ctranspose(wave(:,i)));
        end
    end
    if (timeinterupt=='1')
        % Time evolve through all of the driving cycles
        for z = 1:N
            % Generate the unitary that time evolves the system for each time
            % step of the first driving step
            unitnow = expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/(5*measint)));
            end
            % Iterate over all of the time steps of the first driving step
            for t = 1:measint
                % Time evolve the system for one time step
                density = unitnow*density*ctranspose(unitnow);
                % Iterate over the total number of times that you want entangle
                % a random site with an external qubit as well as the total
                % number of times you want to measure if a particle is located
                % at a particular site.
                for t2i = 1:num
                    % Calculate the eigenvectors V and the eigenvalues I of the
                    % current density matrix
                    [V,I] = eig(density);
                    % Draw a random number
                    draw = rand;
                    % Iterate over all of the entries in probvec
                    for t2 = 1:length(probvec)
                        % If draw is less than the probvec value of the current
                        % iteration, set the y-index value of interest to the
                        % current value of t2.
                        if (draw<probvec(t2))
                            cnow = t2;
                            break;
                        end
                    end
                    % Randomly choose the x-index value of interest
                    ti = randi([0 (Li-1)]);
                    % Randomly choose the value for alpha
                    tk = round(rand)+1;
                    % Calculate the probability for a particle to occupy this
                    % site.
                    probs = abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    % If a random number is less than this probability have the
                    % system only have a population at the site of interest.
                    if (rand<probs)
                        V2 = zeros(2^nqubits);
                        % Iterate over all of the eigenvectors and set the location
                        % corresponding to the site of interest to one.
                        for ti2 = 1:2^nqubits
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = V(tk+2*ti+2*Li*(cnow-1),ti2);
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        density = zeros(2^nqubits);
                        % Reconstruct the density matrix
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    % Otherwise set the system to have a zero probability of
                    % occupying this site and leave the rest of the system
                    % alone
                    else
                        % Store the current eigenvectors
                        V2 = V;
                        % Iterate over all of the eigenvectors
                        for ti2 = 1:2^nqubits
                            % Set the eigenvector of the current iteration to
                            % have a zero probability of occupying the site of
                            % interest.
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = 0;
                            % Normalize this eigenvector.
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    end
                    for t3i = 1:1
                        % Draw a random number
                        Indr = rand;
                        % Iterate over all of the entries in entprob
                        for t2 = 1:length(entprob)
                            % If Indr is less than the current value obtained
                            % from entprob, determine the y-index value of
                            % interest using t2
                            if (Indr<entprob(t2))
                                Indj = t2;
                                break;
                            end
                        end
                        % Randomly choose the x-index value of interest
                        Indi = randi([0 (Li-1)]);
                        % Randomly choose the value for alpha
                        Indk = round(rand)+1;
                        % Add an external qubit
                        density = kron(density,[1 0; 0 0]);
                        % Flip the external qubit if the site of interest is
                        % occupied, otherwise leave it alone.
                        density = measmats(:,:,Indk+2*Indi+2*Li*(Indj-1))*density*ctranspose(measmats(:,:,Indk+2*Indi+2*Li*(Indj-1)));
                        % Remove the external qubit by calculating the reduced
                        % density matrix
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            % Generate the unitary that time evolves the system for each time
            % step of the second driving step
            unitnow = expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/(5*measint)));
            end
            % Iterate over all of the time steps of the second driving step
            for t = 1:measint
                % Time evolve the system for one time step
                density = unitnow*density*ctranspose(unitnow);
                % Iterate over the total number of times that you want entangle
                % a random site with an external qubit as well as the total
                % number of times you want to measure if a particle is located
                % at a particular site.
                for t2i = 1:num
                    % Calculate the eigenvectors V and the eigenvalues I of the
                    % current density matrix
                    [V,I] = eig(density);
                    % Draw a random number
                    draw = rand;
                    % Iterate over all of the entries in probvec
                    for t2 = 1:length(probvec)
                        % If draw is less than the probvec value of the current
                        % iteration, set the y-index value of interest to the
                        % current value of t2.
                        if (draw<probvec(t2))
                            cnow = t2;
                            break;
                        end
                    end
                    % Randomly choose the x-index value of interest
                    ti = randi([0 (Li-1)]);
                    % Randomly choose the value for alpha
                    tk = round(rand)+1;
                    % Calculate the probability for a particle to occupy this
                    % site.
                    probs = abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    % If a random number is less than this probability have the
                    % system only have a population at the site of interest.
                    if (rand<probs)
                        V2 = zeros(2^nqubits);
                        % Iterate over all of the eigenvectors and set the location
                        % corresponding to the site of interest to one.
                        for ti2 = 1:2^nqubits
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = V(tk+2*ti+2*Li*(cnow-1),ti2);
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        density = zeros(2^nqubits);
                        % Reconstruct the density matrix
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    % Otherwise set the system to have a zero probability of
                    % occupying this site and leave the rest of the system
                    % alone
                    else
                        % Store the current eigenvectors
                        V2 = V;
                        % Iterate over all of the eigenvectors
                        for ti2 = 1:2^nqubits
                            % Set the eigenvector of the current iteration to
                            % have a zero probability of occupying the site of
                            % interest.
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = 0;
                            % Normalize this eigenvector.
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    end
                    for t3i = 1:1
                        % Draw a random number
                        Indr = rand;
                        % Iterate over all of the entries in entprob
                        for t2 = 1:length(entprob)
                            % If Indr is less than the current value obtained
                            % from entprob, determine the y-index value of
                            % interest using t2
                            if (Indr<entprob(t2))
                                Indj = t2;
                                break;
                            end
                        end
                        % Randomly choose the x-index value of interest
                        Indi = randi([0 (Li-1)]);
                        % Randomly choose the value for alpha
                        Indk = round(rand)+1;
                        % Add an external qubit
                        density = kron(density,[1 0; 0 0]);
                        % Flip the external qubit if the site of interest is
                        % occupied, otherwise leave it alone.
                        density = measmats(:,:,Indk+2*Indi+2*Li*(Indj-1))*density*ctranspose(measmats(:,:,Indk+2*Indi+2*Li*(Indj-1)));
                        % Remove the external qubit by calculating the reduced
                        % density matrix
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            % Generate the unitary that time evolves the system for each time
            % step of the third driving step
            unitnow = expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/(5*measint)));
            end
            % Iterate over all of the time steps of the third driving step
            for t = 1:measint
                % Time evolve the system for one time step
                density = unitnow*density*ctranspose(unitnow);
                % Iterate over the total number of times that you want entangle
                % a random site with an external qubit as well as the total
                % number of times you want to measure if a particle is located
                % at a particular site.
                for t2i = 1:num
                    % Calculate the eigenvectors V and the eigenvalues I of the
                    % current density matrix
                    [V,I] = eig(density);
                    % Draw a random number
                    draw = rand;
                    % Iterate over all of the entries in probvec
                    for t2 = 1:length(probvec)
                        % If draw is less than the probvec value of the current
                        % iteration, set the y-index value of interest to the
                        % current value of t2.
                        if (draw<probvec(t2))
                            cnow = t2;
                            break;
                        end
                    end
                    % Randomly choose the x-index value of interest
                    ti = randi([0 (Li-1)]);
                    % Randomly choose the value for alpha
                    tk = round(rand)+1;
                    % Calculate the probability for a particle to occupy this
                    % site.
                    probs = abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    % If a random number is less than this probability have the
                    % system only have a population at the site of interest.
                    if (rand<probs)
                        V2 = zeros(2^nqubits);
                        % Iterate over all of the eigenvectors and set the location
                        % corresponding to the site of interest to one.
                        for ti2 = 1:2^nqubits
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = V(tk+2*ti+2*Li*(cnow-1),ti2);
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    % Otherwise set the system to have a zero probability of
                    % occupying this site and leave the rest of the system
                    % alone
                    else
                        % Store the current eigenvectors
                        V2 = V;
                        % Iterate over all of the eigenvectors
                        for ti2 = 1:2^nqubits
                            % Set the eigenvector of the current iteration to
                            % have a zero probability of occupying the site of
                            % interest.
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = 0;
                            % Normalize this eigenvector.
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    end
                    for t3i = 1:1
                        % Draw a random number
                        Indr = rand;
                        % Iterate over all of the entries in entprob
                        for t2 = 1:length(entprob)
                            % If Indr is less than the current value obtained
                            % from entprob, determine the y-index value of
                            % interest using t2
                            if (Indr<entprob(t2))
                                Indj = t2;
                                break;
                            end
                        end
                        % Randomly choose the x-index value of interest
                        Indi = randi([0 (Li-1)]);
                        % Randomly choose the value for alpha
                        Indk = round(rand)+1;
                        % Add an external qubit
                        density = kron(density,[1 0; 0 0]);
                        % Flip the external qubit if the site of interest is
                        % occupied, otherwise leave it alone.
                        density = measmats(:,:,Indk+2*Indi+2*Li*(Indj-1))*density*ctranspose(measmats(:,:,Indk+2*Indi+2*Li*(Indj-1)));
                        % Remove the external qubit by calculating the reduced
                        % density matrix
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            % Generate the unitary that time evolves the system for each time
            % step of the fourth driving step
            unitnow = expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/(5*measint)));
            end
            % Iterate over all of the time steps of the fourth driving step
            for t = 1:measint
                % Time evolve the system for one time step
                density = unitnow*density*ctranspose(unitnow);
                % Iterate over the total number of times that you want entangle
                % a random site with an external qubit as well as the total
                % number of times you want to measure if a particle is located
                % at a particular site.
                for t2i = 1:num
                    % Calculate the eigenvectors V and the eigenvalues I of the
                    % current density matrix
                    [V,I] = eig(density);
                    % Draw a random number
                    draw = rand;
                    % Iterate over all of the entries in probvec
                    for t2 = 1:length(probvec)
                        % If draw is less than the probvec value of the current
                        % iteration, set the y-index value of interest to the
                        % current value of t2.
                        if (draw<probvec(t2))
                            cnow = t2;
                            break;
                        end
                    end
                    % Randomly choose the x-index value of interest
                    ti = randi([0 (Li-1)]);
                    % Randomly choose the value for alpha
                    tk = round(rand)+1;
                    % Calculate the probability for a particle to occupy this
                    % site.
                    probs = abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    % If a random number is less than this probability have the
                    % system only have a population at the site of interest.
                    if (rand<probs)
                        V2 = zeros(2^nqubits);
                        % Iterate over all of the eigenvectors and set the location
                        % corresponding to the site of interest to one.
                        for ti2 = 1:2^nqubits
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = V(tk+2*ti+2*Li*(cnow-1),ti2);
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    % Otherwise set the system to have a zero probability of
                    % occupying this site and leave the rest of the system
                    % alone
                    else
                        % Store the current eigenvectors
                        V2 = V;
                        % Iterate over all of the eigenvectors
                        for ti2 = 1:2^nqubits
                            % Set the eigenvector of the current iteration to
                            % have a zero probability of occupying the site of
                            % interest.
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = 0;
                            % Normalize this eigenvector.
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    end
                    for t3i = 1:1
                        % Draw a random number
                        Indr = rand;
                        % Iterate over all of the entries in entprob
                        for t2 = 1:length(entprob)
                            % If Indr is less than the current value obtained
                            % from entprob, determine the y-index value of
                            % interest using t2
                            if (Indr<entprob(t2))
                                Indj = t2;
                                break;
                            end
                        end
                        % Randomly choose the x-index value of interest
                        Indi = randi([0 (Li-1)]);
                        % Randomly choose the value for alpha
                        Indk = round(rand)+1;
                        % Add an external qubit
                        density = kron(density,[1 0; 0 0]);
                        % Flip the external qubit if the site of interest is
                        % occupied, otherwise leave it alone.
                        density = measmats(:,:,Indk+2*Indi+2*Li*(Indj-1))*density*ctranspose(measmats(:,:,Indk+2*Indi+2*Li*(Indj-1)));
                        % Remove the external qubit by calculating the reduced
                        % density matrix
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            %%%
            % Generate the unitary that time evolves the system for each time
            % step of the fifth driving step
            unitnow = expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/(5*measint));
            for t = 2:ntimes
                unitnow = kron(unitnow,expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/(5*measint)));
            end
            % Iterate over all of the time steps of the fifth driving step
            for t = 1:measint
                % Time evolve the system for one time step
                density = unitnow*density*ctranspose(unitnow);
                % Iterate over the total number of times that you want entangle
                % a random site with an external qubit as well as the total
                % number of times you want to measure if a particle is located
                % at a particular site.
                for t2i = 1:num
                    % Calculate the eigenvectors V and the eigenvalues I of the
                    % current density matrix
                    [V,I] = eig(density);
                    % Draw a random number
                    draw = rand;
                    % Iterate over all of the entries in probvec
                    for t2 = 1:length(probvec)
                        % If draw is less than the probvec value of the current
                        % iteration, set the y-index value of interest to the
                        % current value of t2.
                        if (draw<probvec(t2))
                            cnow = t2;
                            break;
                        end
                    end
                    % Randomly choose the x-index value of interest
                    ti = randi([0 (Li-1)]);
                    % Randomly choose the value for alpha
                    tk = round(rand)+1;
                    % Calculate the probability for a particle to occupy this
                    % site.
                    probs = abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    % If a random number is less than this probability have the
                    % system only have a population at the site of interest.
                    if (rand<probs)
                        V2 = zeros(2^nqubits);
                        % Iterate over all of the eigenvectors and set the location
                        % corresponding to the site of interest to one.
                        for ti2 = 1:2^nqubits
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = V(tk+2*ti+2*Li*(cnow-1),ti2);
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    % Otherwise set the system to have a zero probability of
                    % occupying this site and leave the rest of the system
                    % alone
                    else
                        % Store the current eigenvectors
                        V2 = V;
                        % Iterate over all of the eigenvectors
                        for ti2 = 1:2^nqubits
                            % Set the eigenvector of the current iteration to
                            % have a zero probability of occupying the site of
                            % interest.
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = 0;
                            % Normalize this eigenvector.
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    end
                    for t3i = 1:1
                        % Draw a random number
                        Indr = rand;
                        % Iterate over all of the entries in entprob
                        for t2 = 1:length(entprob)
                            % If Indr is less than the current value obtained
                            % from entprob, determine the y-index value of
                            % interest using t2
                            if (Indr<entprob(t2))
                                Indj = t2;
                                break;
                            end
                        end
                        % Randomly choose the x-index value of interest
                        Indi = randi([0 (Li-1)]);
                        % Randomly choose the value for alpha
                        Indk = round(rand)+1;
                        % Add an external qubit
                        density = kron(density,[1 0; 0 0]);
                        % Flip the external qubit if the site of interest is
                        % occupied, otherwise leave it alone.
                        density = measmats(:,:,Indk+2*Indi+2*Li*(Indj-1))*density*ctranspose(measmats(:,:,Indk+2*Indi+2*Li*(Indj-1)));
                        % Remove the external qubit by calculating the reduced
                        % density matrix
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
            end
            % Calculate the probability for the particle to occupy each of the
            % sites
            for j = 1:num
                jprobsb(1,j,z) = abs(density(j,j));
            end
        end
    else
        % Calculate after how many driving steps, the entanglement and wave
        % function collapse occurs
        measint2 = round(1/measint);
        aph = 0;
        % Iterate over all driving cycles
        for z = 1:N
            % Iterate over all driving steps
            for z2 = 1:5
                aph = aph + 1;
                % Implement the first driving step if z2==1
                if (z2==1)
                    unitnow = expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % Implement the second driving step if z2==2
                elseif (z2==2)
                    unitnow = expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % Implement the third driving step if z2==3
                elseif (z2==3)
                    unitnow = expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % Implement the fourth driving step if z2==4
                elseif (z2==4)
                    unitnow = expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                % Implement the fifth driving step if z2==5
                elseif (z2==5)
                    unitnow = expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/5);
                    for z3 = 2:ntimes
                        unitnow = kron(unitnow,expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/5));
                    end
                    density = unitnow*density*ctranspose(unitnow);
                end
                % After the appropriate driving steps, implement the
                % entanglement and wave function collapse methods
                if (mod(aph,measint2)==0)
                % Iterate over the number of times that we want to implement
                % the entanglement and wave function collapse operations
                for t2i = 1:num
                    % Calculate the eigenvectors V and the eigenvalues I of the
                    % current density matrix
                    [V,I] = eig(density);
                    % Draw a random number
                    draw = rand;
                    % Iterate over all of the entries in probvec
                    for t2 = 1:length(probvec)
                        % If draw is less than the probvec value of the current
                        % iteration, set the y-index value of interest to the
                        % current value of t2.
                        if (draw<probvec(t2))
                            cnow = t2;
                            break;
                        end
                    end
                    % Randomly choose the x-index value of interest
                    ti = randi([0 (Li-1)]);
                    % Randomly choose the value for alpha
                    tk = round(rand)+1;
                    % Calculate the probability for a particle to occupy this
                    % site.
                    probs = abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    % If a random number is less than this probability have the
                    % system only have a population at the site of interest.
                    if (rand<probs)
                        V2 = zeros(2^nqubits);
                        % Iterate over all of the eigenvectors and set the location
                        % corresponding to the site of interest to one.
                        for ti2 = 1:2^nqubits
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = V(tk+2*ti+2*Li*(cnow-1),ti2);
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    % Otherwise set the system to have a zero probability of
                    % occupying this site and leave the rest of the system
                    % alone
                    else
                        % Store the current eigenvectors
                        V2 = V;
                        % Iterate over all of the eigenvectors
                        for ti2 = 1:2^nqubits
                            % Set the eigenvector of the current iteration to
                            % have a zero probability of occupying the site of
                            % interest.
                            V2(tk+2*ti+2*Li*(cnow-1),ti2) = 0;
                            % Normalize this eigenvector.
                            norm = abs(ctranspose(V2(:,ti2))*V2(:,ti2));
                            if (norm>0)
                                V2(:,ti2) = V2(:,ti2)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti2 = 1:2^nqubits
                            density = density + I(ti2,ti2)*V2(:,ti2)*ctranspose(V2(:,ti2));
                        end
                        density = density/trace(abs(density));
                    end
                    for t3i = 1:1
                        % Draw a random number
                        Indr = rand;
                        % Iterate over all of the entries in entprob
                        for t2 = 1:length(entprob)
                            % If Indr is less than the current value obtained
                            % from entprob, determine the y-index value of
                            % interest using t2
                            if (Indr<entprob(t2))
                                Indj = t2;
                                break;
                            end
                        end
                        % Randomly choose the x-index value of interest
                        Indi = randi([0 (Li-1)]);
                        % Randomly choose the value for alpha
                        Indk = round(rand)+1;
                        % Add an external qubit
                        density = kron(density,[1 0; 0 0]);
                        % Flip the external qubit if the site of interest is
                        % occupied, otherwise leave it alone.
                        density = measmats(:,:,Indk+2*Indi+2*Li*(Indj-1))*density*ctranspose(measmats(:,:,Indk+2*Indi+2*Li*(Indj-1)));
                        % Remove the external qubit by calculating the reduced
                        % density matrix
                        [rdensity] = ReducedDensity(density,ntimes*nqubits+1,1:(ntimes*nqubits));
                        density = rdensity;
                    end
                end
                end
                % Calculate the probability for the particle to occupy each of
                % the sites
                if (z2==5)
                    for j = 1:num
                        jprobsb(1,j,z) = abs(density(j,j));
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
