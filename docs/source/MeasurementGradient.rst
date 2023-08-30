=========================================================================================================================================
Presentation of the Code that Tests the Model that Uses a Gradient that Determines How Often the Particle has its Wave Function Collapsed
=========================================================================================================================================

The code presented below tests the model where the likelihood of the particle located at a particular site having its wave function collapse is dependent on the value of the site's y-index. The higher the value of the y-index (starting at zero), the more likely the particle is to have its wave function collapse. This works by dividing the time evolution of any particular driving step into a series of time steps. After each time step, a y-index is randomly chosen where y-indices that have a higher value are more likely to be chosen than ones that have a lower y-index. Then the eigenvectors and eigenvalues of the relevant density matrix are calculated as well as the probability for the particle to occupy the chosen y-index. A random number generator then decides using the probability calculated if all of the components of the eigenvectors corresponding to the chosen y-index are set to zero or if all of the other components are set to zero. The eigenvectors are then normalized and a new density matrix is reconstructed using the eigenvalues and the new eigenvectors. Finally, the density matrix itself is normalized. The main MATLAB file that runs the sub-files that obtain the data and plots the resulting data is given below:

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
   % Store the probability of the particle occupying each of the y-indices for
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
   % Calculate the standard error associated with the particle occupying
   % each of the y-indices for each of the driving cycles for the system
   % unaffected by wave function collapse operations
   jprobs1std = zeros(NTime,Lj);
   % Calculate the standard error associated with the particle occupying
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
    % Determine the number of times per driving step that the presence of a
    % particle is measured for a single site
    measint = 100;
    % The following if else statements determines how the time evolution takes
    % place depending on the size of measint
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
    % steps. The velocity matrices V1 and V3 are irrelevent for this
    % particular calculation.
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
    % particle occupying each site for the reference system unaffected by wave
    % function collapse
    jprobsa = zeros(1,2+2*(Li-1)+2*Li*(Lj-1),N);
    % The following matrix stores information regarding the probability of the
    % particle occupying each site for the system where wave function collapse
    % is involved
    jprobsb = zeros(1,2+2*(Li-1)+2*Li*(Lj-1),N);
    aph = 0;
    % The following matrix stores the projection operators that are used to
    % calculate the probability of the particle occupying each of the sites
    sitexpectations = zeros(2^(ntimes*nqubits),2^(ntimes*nqubits),2+2*(Li-1)+2*Li*(Lj-1));
    % The following pretty much does the same thing but adds an extra dimension
    % to sort the projection operators according to the y-indices. jvales is actually
    % used to store the wave functions classified by the relevant y-index.
    jvals = zeros(2^(ntimes*nqubits),round(2^(ntimes*nqubits)/Lj),Lj);
    for j = 0:(Lj-1)
        saph = 0;
        for i = 0:(Li-1)
            for k = 1:2
                aph = aph + 1;
                saph = saph + 1;
                sitexpectations(k+2*i+2*Li*j,k+2*i+2*Li*j,aph) = 1;
                jvals(k+2*i+2*Li*j,saph,j+1) = 1;
            end
        end
    end
    % Stores how many sites are in the system
    num = aph;
    % Time evolve the system that is unaffected by wave function collapse
    % operations
    for z = 1:N
        wave2 = expm(-1i*(H5)*(1+TimeDisorder5(z))*2*pi/5)*expm(-1i*(H4)*(1+TimeDisorder4(z))*2*pi/5)*expm(-1i*(H3)*(1+TimeDisorder3(z))*2*pi/5)*expm(-1i*(H2)*(1+TimeDisorder2(z))*2*pi/5)*expm(-1i*(H1)*(1+TimeDisorder1(z))*2*pi/5)*wave2;
        % Calculate the probability for the particle to occupy each of the
        % sites
        for j = 1:num
            jprobsa(1,j,z) = ctranspose(wave2)*sitexpectations(:,:,j)*wave2;
        end
    end
    % Generate the density matrix for the system where wave function collapse
    % operations are involved.
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
                % Calculate the eigenvectors V and the eigenvalues I of the
                % current density matrix
                [V,I] = eig(density);
                % Draw a random number
                draw = rand;
                % Iterate over all of the entries in probvec
                for t2 = 1:length(probvec)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<probvec(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Store the probability of the particle to occupy a given value
                % of the y-index
                probs = 0;
                % Iterate over all values of the x-index for the given y-index
                for ti = 0:(Li-1)
                    % Iterate over all values of alpha
                    for tk = 1:2
                        probs = probs + abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    end
                end
                % If a random number is less than the probability for the
                % particle to occupy a y-index of interest
                if (rand<probs)
                    V2 = zeros(2^nqubits);
                    % Iterate over all of the eigenvectors V
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites that correspond to the
                        % y-index of interest
                        for tj = 1:round((2^nqubits)/Lj)
                            % Populate V2 with entries that preserve the
                            % information of the sites with the y-index of
                            % interest
                            V2(:,ti) = V2(:,ti) + ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the corresponding vectors
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
                else
                    % Copy the information of the eigenvectors
                    V2 = V;
                    % Iterate over all of the eigenvectors
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites for any given value of
                        % alpha or the x-index
                        for tj = 1:round((2^nqubits)/Lj)
                            % Erase the information of the site that is
                            % currently being iterated over, such that it is
                            % now zero
                            V2(:,ti) = V2(:,ti) - ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the resulting vector
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
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
                % Calculate the eigenvectors V and the eigenvalues I of the
                % current density matrix
                [V,I] = eig(density);
                % Draw a random number
                draw = rand;
                % Iterate over all of the entries in probvec
                for t2 = 1:length(probvec)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<probvec(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Store the probability of the particle to occupy a given value
                % of the y-index
                probs = 0;
                % Iterate over all values of the x-index for the given y-index
                for ti = 0:(Li-1)
                    % Iterate over all values of alpha
                    for tk = 1:2
                        probs = probs + abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    end
                end
                % If a random number is less than the probability for the
                % particle to occupy a y-index of interest
                if (rand<probs)
                    V2 = zeros(2^nqubits);
                    % Iterate over all of the eigenvectors V
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites that correspond to the
                        % y-index of interest
                        for tj = 1:round((2^nqubits)/Lj)
                            % Populate V2 with entries that preserve the
                            % information of the sites with the y-index of
                            % interest
                            V2(:,ti) = V2(:,ti) + ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the corresponding vectors
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
                else
                    % Copy the information of the eigenvectors
                    V2 = V;
                    % Iterate over all of the eigenvectors
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites for any given value of
                        % alpha or the x-index
                        for tj = 1:round((2^nqubits)/Lj)
                            % Erase the information of the site that is
                            % currently being iterated over, such that it is
                            % now zero
                            V2(:,ti) = V2(:,ti) - ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the resulting vector
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
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
                % Calculate the eigenvectors V and the eigenvalues I of the
                % current density matrix
                [V,I] = eig(density);
                % Draw a random number
                draw = rand;
                % Iterate over all of the entries in probvec
                for t2 = 1:length(probvec)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<probvec(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Store the probability of the particle to occupy a given value
                % of the y-index
                probs = 0;
                % Iterate over all values of the x-index for the given y-index
                for ti = 0:(Li-1)
                    % Iterate over all values of alpha
                    for tk = 1:2
                        probs = probs + abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    end
                end
                % If a random number is less than the probability for the
                % particle to occupy a y-index of interest
                if (rand<probs)
                    V2 = zeros(2^nqubits);
                    % Iterate over all of the eigenvectors V
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites that correspond to the
                        % y-index of interest
                        for tj = 1:round((2^nqubits)/Lj)
                            % Populate V2 with entries that preserve the
                            % information of the sites with the y-index of
                            % interest
                            V2(:,ti) = V2(:,ti) + ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the corresponding vectors
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
                else
                    % Copy the information of the eigenvectors
                    V2 = V;
                    % Iterate over all of the eigenvectors
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites for any given value of
                        % alpha or the x-index
                        for tj = 1:round((2^nqubits)/Lj)
                            % Erase the information of the site that is
                            % currently being iterated over, such that it is
                            % now zero
                            V2(:,ti) = V2(:,ti) - ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the resulting vector
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
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
                % Calculate the eigenvectors V and the eigenvalues I of the
                % current density matrix
                [V,I] = eig(density);
                % Draw a random number
                draw = rand;
                % Iterate over all of the entries in probvec
                for t2 = 1:length(probvec)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<probvec(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Store the probability of the particle to occupy a given value
                % of the y-index
                probs = 0;
                % Iterate over all values of the x-index for the given y-index
                for ti = 0:(Li-1)
                    % Iterate over all values of alpha
                    for tk = 1:2
                        probs = probs + abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    end
                end
                % If a random number is less than the probability for the
                % particle to occupy a y-index of interest
                if (rand<probs)
                    V2 = zeros(2^nqubits);
                    % Iterate over all of the eigenvectors V
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites that correspond to the
                        % y-index of interest
                        for tj = 1:round((2^nqubits)/Lj)
                            % Populate V2 with entries that preserve the
                            % information of the sites with the y-index of
                            % interest
                            V2(:,ti) = V2(:,ti) + ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the corresponding vectors
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
                else
                    % Copy the information of the eigenvectors
                    V2 = V;
                    % Iterate over all of the eigenvectors
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites for any given value of
                        % alpha or the x-index
                        for tj = 1:round((2^nqubits)/Lj)
                            % Erase the information of the site that is
                            % currently being iterated over, such that it is
                            % now zero
                            V2(:,ti) = V2(:,ti) - ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the resulting vector
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
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
                % Calculate the eigenvectors V and the eigenvalues I of the
                % current density matrix
                [V,I] = eig(density);
                % Draw a random number
                draw = rand;
                % Iterate over all of the entries in probvec
                for t2 = 1:length(probvec)
                    % If draw is less than the probvec value of the current
                    % iteration, set the y-index value of interest according to
                    % the current value of t2.
                    if (draw<probvec(t2))
                        cnow = t2;
                        break;
                    end
                end
                % Store the probability of the particle to occupy a given value
                % of the y-index
                probs = 0;
                % Iterate over all values of the x-index for the given y-index
                for ti = 0:(Li-1)
                    % Iterate over all values of alpha
                    for tk = 1:2
                        probs = probs + abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                    end
                end
                % If a random number is less than the probability for the
                % particle to occupy a y-index of interest
                if (rand<probs)
                    V2 = zeros(2^nqubits);
                    % Iterate over all of the eigenvectors V
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites that correspond to the
                        % y-index of interest
                        for tj = 1:round((2^nqubits)/Lj)
                            % Populate V2 with entries that preserve the
                            % information of the sites with the y-index of
                            % interest
                            V2(:,ti) = V2(:,ti) + ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the corresponding vectors
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
                else
                    % Copy the information of the eigenvectors
                    V2 = V;
                    % Iterate over all of the eigenvectors
                    for ti = 1:2^nqubits
                        % Iterate over all of the sites for any given value of
                        % alpha or the x-index
                        for tj = 1:round((2^nqubits)/Lj)
                            % Erase the information of the site that is
                            % currently being iterated over, such that it is
                            % now zero
                            V2(:,ti) = V2(:,ti) - ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                        end
                        % Normalize the resulting vector
                        norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                        if (norm>0)
                            V2(:,ti) = V2(:,ti)/sqrt(norm);
                        end
                    end
                    % Reconstruct the density matrix
                    density = zeros(2^nqubits);
                    for ti = 1:2^nqubits
                        density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                    end
                    density = density/trace(abs(density));
                end
            end
            % Calculate the probability for the particle to occupy each of the
            % sites
            for j = 1:num
                jprobsb(1,j,z) = abs(density(j,j));
            end     
        end
    else
        % Calculate after how many driving steps the wave function collapse
        % operations are implemented
        measint2 = round(1/measint);
        aph = 0;
        % Iterate over all of the driving cycles
        for z = 1:N
            % Iterate over all of the driving steps
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
                % After the appropriate number of driving steps, implement the
                % wave function collapse operations
                if (mod(aph,measint2)==0)
                    % Calculate the eigenvectors V and the eigenvalues I of the
                    % current density matrix
                    [V,I] = eig(density);
                    % Draw a random number
                    draw = rand;
                    % Iterate over all of the entries in probvec
                    for t2 = 1:length(probvec)
                        % If draw is less than the probvec value of the current
                        % iteration, set the y-index value of interest
                        % according to the current value of t2.
                        if (draw<probvec(t2))
                            cnow = t2;
                            break;
                        end
                    end
                    % Store the probability of the particle to occupy a given
                    % value of the y-index
                    probs = 0;
                    % Iterate over all values of the x-index for the given
                    % y-index
                    for ti = 0:(Li-1)
                        % Iterate over all values of alpha
                        for tk = 1:2
                            probs = probs + abs(density(tk+2*ti+2*Li*(cnow-1),tk+2*ti+2*Li*(cnow-1)));
                        end
                    end
                    % If a random number is less than the probability for the
                    % particle to occupy a y-index of interest
                    if (rand<probs)
                        V2 = zeros(2^nqubits);
                        % Iterate over all of the eigenvectors V
                        for ti = 1:2^nqubits
                            % Iterate over all of the sites that correspond to
                            % the y-index of interest
                            for tj = 1:round((2^nqubits)/Lj)
                                % Populate V2 with entries that preserve the
                                % information of the sites with the y-index of
                                % interest
                                V2(:,ti) = V2(:,ti) + ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                            end
                            % Normalize the corresponding vectors
                            norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                            if (norm>0)
                                V2(:,ti) = V2(:,ti)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti = 1:2^nqubits
                            density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                        end
                        density = density/trace(abs(density));
                    else
                        % Copy the information of the eigenvectors
                        V2 = V;
                        % Iterate over all of the eigenvectors
                        for ti = 1:2^nqubits
                            % Iterate over all of the sites for any given value of
                            % alpha or the x-index
                            for tj = 1:round((2^nqubits)/Lj)
                                % Erase the information of the site that is
                                % currently being iterated over, such that it
                                % is now zero
                                V2(:,ti) = V2(:,ti) - ctranspose(jvals(:,tj,cnow))*V(:,ti)*jvals(:,tj,cnow);
                            end
                            % Normalize the resulting vector
                            norm = abs(ctranspose(V2(:,ti))*V2(:,ti));
                            if (norm>0)
                                V2(:,ti) = V2(:,ti)/sqrt(norm);
                            end
                        end
                        % Reconstruct the density matrix
                        density = zeros(2^nqubits);
                        for ti = 1:2^nqubits
                            density = density + I(ti,ti)*V2(:,ti)*ctranspose(V2(:,ti));
                        end
                        density = density/trace(abs(density));
                    end
                end
                % At the end of each driving step, calculate the probability
                % for the particle to occupy each of the sites.
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
    % appropriate locations such that they perform the actions they were
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
