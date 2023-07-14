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

The code above runs TwoDimxyQ.m, which is the main file that actually runs the simulation for each noise realization. This code is presented below:

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
