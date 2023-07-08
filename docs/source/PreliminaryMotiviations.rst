======================================================
Preliminary Motivations for 'Quantum Gravity' Research
======================================================

This documentation prevents, as the title says, the preliminary motivations for this research on quantum gravity. This code allows for the storage of a set of single qubits in a vector that is much more compressed than it otherwise would be. This works by having a series of control qubits and a single target qubit that store the information of these single qubits. If there is a single control qubit, then the total wave function can store the information of two single qubits according to:

:math:`$$\Psi_{\mathrm{tot}}^T = [1\ 0]\otimes [a\ b] + [0\ 1]\otimes [c\ d] = [a\ b\ c\ d]$$`

Meanwhile, if there are two control qubits, the total wave function becomes:

:math:`$$\Psi_{\mathrm{tot}}^T = [1\ 0\ 0\ 0]\otimes [a\ b] + [0\ 1\ 0\ 0]\otimes [c\ d] + [0\ 0\ 1\ 0]\otimes [e\ f] + [0\ 0\ 0\ 1]\otimes [g\ h] = [a\ b\ c\ d\ e\ f\ g\ h]$$`

Therefore, as the number of 'compressed' qubits (n) stored in this wave function grows by a factor of two, the size of wave function grows by a factor of :math:`$2\times 2\times n$` as opposed to the normal :math:`$2^{n}$` for the normal storage of qubits in a wave function. The other interesting aspect to this is that as the number control qubits grows, the number of shots needed to gain the same accuracy of the measurement results grows by a factor two. So while, in a normal ket of unentangled qubits where the number of shots needed to obtain a certain degree of accuracy scales by a factor of one (whatever that means), the number of shots needed to obtain the same accuracy with these 'compressed' qubits scales as a factor of :math:`$2^n$`. This can be seen as a kind of length contraction, where the wave function shortens to becoming :math:`$2\times 2\times n$` as opposed to :math:`$2^{n}$`, as well as a kind of time dilation, where the number of shots needed for a certain degree of accuracy increases by a factor of :math:`$2^n$` as opposed to the 'uncompressed qubits' that only increase by a factor of 'one'. The code for this is presented below:

.. code-block:: matlab

   clear;
   clc;
   % 2 to the power nqubits determines how many 'compressed' and how many
   % 'uncompressed' qubits we are comparing
   nqubits = 4;
   % Store all of the eigenvectors of the Pauli matrices in spinvecs
   spinvecs = zeros(2,2,3);
   spinvecs(:,1,1) = [1; 0];
   spinvecs(:,2,1) = [0; 1];
   spinvecs(:,1,2) = [1; -1]/sqrt(2);
   spinvecs(:,2,2) = [1; 1]/sqrt(2);
   spinvecs(:,1,3) = [1; -1i]/sqrt(2);
   spinvecs(:,2,3) = [1; 1i]/sqrt(2);
   rng('shuffle')
   % Define how both the compressed and uncompressed qubits will rotate around
   % the x-axis
   thetaxs = 2*pi*rand(1,2^nqubits);
   % Define how both the compressed and uncompressed qubits will rotate around
   % the z-axis
   thetazs = 2*pi*rand(1,2^nqubits);
   % Generate the rotation matrices about the x and z axes
   rotationxs = zeros(2,2,2^nqubits);
   rotationzs = zeros(2,2,2^nqubits);
   for i = 1:2^nqubits
       rotationxs(:,:,i) = expm(-1i*thetaxs(i)*[0 1; 1 0]/2);
       rotationzs(:,:,i) = expm(-1i*thetazs(i)*[1 0; 0 -1]/2);
   end
   % Rotate the uncompressed qubits about the x and z axes and store them
   % within qubits1
   qubits1 = zeros(2,2,2^nqubits);
   for i = 1:2^nqubits
       qubits1(:,:,i) = rotationzs(:,:,i)*rotationxs(:,:,i)*[1 0; 0 0]*ctranspose(rotationxs(:,:,i))*ctranspose(rotationzs(:,:,i));
   end
   % Reorient the compressed qubits such that it becomes very difficult for an
   % outside party to guess the orientation
   key = randperm(2^nqubits);
   % Randomly determine whether the control qubits that are used to store the
   % compressed qubits are oriented in the x, y, or z directions
   basis = randi([1 3],1,nqubits);
   % If the relevant basis is oriented in the x-direction, have a Hadamard
   % gate act on the qubit, otherwise leave the qubit alone
   if (basis(1)==1)
       initrot = [1 1; 1 -1]/sqrt(2);
   else
       initrot = eye(2);
   end
   for i = 2:nqubits
       if (basis(i)==1)
           initrot = kron(initrot,[1 1; 1 -1]/sqrt(2));
       else
           initrot = kron(initrot,eye(2));
       end
   end
   initrot = kron(initrot,eye(2));
   % Set the initial state for the compressed qubits
   qubits2 = [1 0; 0 0];
   for i = 1:nqubits
       qubits2 = kron(qubits2,[1 0; 0 0]);
   end
   % Rotate the control qubits, which are used to store the compressed qubits,
   % into the appropriate basis before the actual control operations are
   % performed
   qubits2 = initrot*qubits2*ctranspose(initrot);
   % Construct the transformation that encodes the compressed qubits onto the
   % target qubit
   transform = zeros(2^(nqubits+1));
   % Iterate through all of the compressed qubits
   for i = 1:2^nqubits
       % Turn the number of the current iteration into a bitstring of the
       % appropriate length
       const = dec2bin(i-1);
       const2 = nqubits - length(const);
       for j = 1:const2
           const = ['0' const];
       end
       % Construct the control operations in the appropriate basis for
       % encryption purposes
       if (const(1)=='1')
           rotation = spinvecs(:,2,basis(1))*ctranspose(spinvecs(:,2,basis(1)));
       else
           rotation = spinvecs(:,1,basis(1))*ctranspose(spinvecs(:,1,basis(1)));
       end
       for j = 2:nqubits
           if (const(j)=='1')
               rotation = kron(rotation,spinvecs(:,2,basis(j))*ctranspose(spinvecs(:,2,basis(j))));
           else
               rotation = kron(rotation,spinvecs(:,1,basis(j))*ctranspose(spinvecs(:,1,basis(j))));
           end
       end
       % Encode the compressed qubit considering the appropriate orientation
       % of the control qubit as given by the key
       rotation = kron(rotation,rotationzs(:,:,key(i))*rotationxs(:,:,key(i)));
       transform = transform + rotation;
   end
   % Have the control operations act on the sets of qubits that store the
   % compressed qubits
   qubits2 = transform*qubits2*ctranspose(transform);
   % Determine the number of shots for the uncompressed qubits
   shots1 = 100000;
   % The number of shots for the compressed qubits is set to be 2 to the power
   % of nqubits multiplied by the number of shots for the uncompressed qubits
   % in order to obtain a similar level of accuracy
   shots2 = 2^(nqubits)*shots1;
   % Determine the probability of the sampled outcomes for measurements in the
   % z-direction for the uncompressed qubits
   probabilityz1 = zeros(1,2^nqubits);
   % Determine the probability of the sampled outcomes for measurements in the
   % x-direction for the uncompressed qubits
   probabilityx1 = zeros(1,2^nqubits);
   % Determine the probability of the sampled outcomes for measurements in the
   % y-direction for the uncompressed qubits
   probabilityy1 = zeros(1,2^nqubits);
   for i = 1:2^nqubits
       probabilityz1(i) = abs(trace(qubits1(:,:,i)*spinvecs(:,1,1)*ctranspose(spinvecs(:,1,1))));
       probabilityx1(i) = abs(trace(qubits1(:,:,i)*spinvecs(:,1,2)*ctranspose(spinvecs(:,1,2))));
       probabilityy1(i) = abs(trace(qubits1(:,:,i)*spinvecs(:,1,3)*ctranspose(spinvecs(:,1,3))));
   end
   % Determine the probability of the sampled outcomes for measurements in the
   % z-direction for the compressed qubits
   probabilityz2 = zeros(2,2^nqubits);
   % Determine the probability of the sampled outcomes for measurements in the
   % x-direction for the compressed qubits
   probabilityx2 = zeros(2,2^nqubits);
   % Determine the probability of the sampled outcomes for measurements in the
   % y-direction for the compressed qubits
   probabilityy2 = zeros(2,2^nqubits);
   probz = 0;
   probx = 0;
   proby = 0;
   for i = 1:2^nqubits
       const = dec2bin(i-1);
       const2 = nqubits - length(const);
       for j = 1:const2
           const = ['0' const];
       end
       % Use the appropriate bases for the control qubits
       if (const(1)=='1')
           rotation = spinvecs(:,2,basis(1))*ctranspose(spinvecs(:,2,basis(1)));
       else
           rotation = spinvecs(:,1,basis(1))*ctranspose(spinvecs(:,1,basis(1)));
       end
       for j = 2:nqubits
           if (const(j)=='1')
               rotation = kron(rotation,spinvecs(:,2,basis(j))*ctranspose(spinvecs(:,2,basis(j))));
           else
               rotation = kron(rotation,spinvecs(:,1,basis(j))*ctranspose(spinvecs(:,1,basis(j))));
           end
       end
       % Iterate over the spin up and down directions for the relevant bases
       for j = 1:2
           probz = probz + abs(trace(qubits2*kron(rotation,spinvecs(:,j,1)*ctranspose(spinvecs(:,j,1)))));
           probabilityz2(j,i) = probz;
           probx = probx + abs(trace(qubits2*kron(rotation,spinvecs(:,j,2)*ctranspose(spinvecs(:,j,2)))));
           probabilityx2(j,i) = probx;
           proby = proby + abs(trace(qubits2*kron(rotation,spinvecs(:,j,3)*ctranspose(spinvecs(:,j,3)))));
           probabilityy2(j,i) = proby;
       end
   end
   % Count the sampled outcomes in the z-direction for the uncompressed qubits
   samplez1 = zeros(2,2^nqubits);
   % Count the sampled outcomes in the x-direction for the uncompressed qubits
   samplex1 = zeros(2,2^nqubits);
   % Count the sampled outcomes in the y-direction for the uncompressed qubits
   sampley1 = zeros(2,2^nqubits);
   rng('shuffle')
   % Iterate over all of the uncompressed qubits
   for i = 1:2^nqubits
       % Iterate over all of the shots
       for j = 1:shots1
           % Have a random number generator randomly choose the outcome in the
           % z-direction
           const = rand;
           if (const<probabilityz1(i))
               samplez1(1,i) = samplez1(1,i) + 1;
           else
               samplez1(2,i) = samplez1(2,i) + 1;
           end
           % Have a random number generator randomly choose the outcome in the
           % x-direction
           const = rand;
           if (const<probabilityx1(i))
               samplex1(1,i) = samplex1(1,i) + 1;
           else
               samplex1(2,i) = samplex1(2,i) + 1;
           end
           % Have a random number generator randomly choose the outcome in the
           % y-direction
           const = rand;
           if (const<probabilityy1(i))
               sampley1(1,i) = sampley1(1,i) + 1;
           else
               sampley1(2,i) = sampley1(2,i) + 1;
           end
       end
   end
   rng('shuffle')
   % Count the sampled outcomes in the z-direction for the compressed qubits
   samplez2i = zeros(2,2^nqubits);
   % Count the sampled outcomes in the x-direction for the compressed qubits
   samplex2i = zeros(2,2^nqubits);
   % Count the sampled outcomes in the y-direction for the compressed qubits
   sampley2i = zeros(2,2^nqubits);
   % Count the number of times that the wave function collapses to a
   % particular configuration for the compressed qubits; where the outcome of
   % the target qubit is measured in the z-direction. This forms an effective
   % shot count.
   shotcountzi = zeros(1,2^nqubits);
   % Count the number of times that the wave function collapses to a
   % particular configuration for the compressed qubits; where the outcome of
   % the target qubit is measured in the x-direction. This forms an effective
   % shot count.
   shotcountxi = zeros(1,2^nqubits);
   % Count the number of times that the wave function collapses to a
   % particular configuration for the compressed qubits; where the outcome of
   % the target qubit is measured in the y-direction. This forms an effective
   % shot count.
   shotcountyi = zeros(1,2^nqubits);
   % Iterate over all of the compressed qubits
   for i = 1:shots2
       % Have a random number generator determine the outcome
       const = rand;
       aph = 0;
       % Iterate over all of the probability distributions to determine where
       % the outcome lands
       for j = 1:2^nqubits
           for k = 1:2
               if (const<probabilityz2(k,j))
                   % Count the number of sampled outcomes in the z-direction
                   samplez2i(k,j) = samplez2i(k,j) + 1;
                   % Count the number of effective shot counts in the z-direction
                   shotcountzi(j) = shotcountzi(j) + 1;
                   aph = 1;
                   break;
               end
           end
           if (aph==1)
               break;
           end
       end
       % Repeat the process for the x-direction
       const = rand;
       aph = 0;
       for j = 1:2^nqubits
           for k = 1:2
               if (const<probabilityx2(k,j))
                   samplex2i(k,j) = samplex2i(k,j) + 1;
                   shotcountxi(j) = shotcountxi(j) + 1;
                   aph = 1;
                   break;
               end
           end
           if (aph==1)
               break;
           end
       end
       % Repeat the process for the y-direction
       const = rand;
       aph = 0;
       for j = 1:2^nqubits
           for k = 1:2
               if (const<probabilityy2(k,j))
                   sampley2i(k,j) = sampley2i(k,j) + 1;
                   shotcountyi(j) = shotcountyi(j) + 1;
                   aph = 1;
                   break;
               end
           end
           if (aph==1)
               break;
           end
       end
   end
   % Reorient the sampled outcomes and the effective shot counts according to
   % the key in order to decrypt the compressed qubits.
   samplez2 = zeros(2,2^nqubits);
   samplex2 = zeros(2,2^nqubits);
   sampley2 = zeros(2,2^nqubits);
   shotcountz = zeros(1,2^nqubits);
   shotcountx = zeros(1,2^nqubits);
   shotcounty = zeros(1,2^nqubits);
   for i = 1:2^nqubits
       samplez2(:,i) = samplez2i(:,find(key==i));
       shotcountz(i) = shotcountzi(find(key==i));
       samplex2(:,i) = samplex2i(:,find(key==i));
       shotcountx(i) = shotcountxi(find(key==i));
       sampley2(:,i) = sampley2i(:,find(key==i));
       shotcounty(i) = shotcountyi(find(key==i));
   end
   probz1 = samplez1/shots1;
   probx1 = samplex1/shots1;
   proby1 = sampley1/shots1;
   probz2 = zeros(2,2^nqubits);
   probx2 = zeros(2,2^nqubits);
   proby2 = zeros(2,2^nqubits);
   for i = 1:2^nqubits
       probz2(:,i) = samplez2(:,i)/shotcountz(i);
       probx2(:,i) = samplex2(:,i)/shotcountx(i);
       proby2(:,i) = sampley2(:,i)/shotcounty(i);
   end
   % Display all of the data
   disp('Sampled probability for the uncompressed qubits in the z-down direction is:')
   disp(probz1(1,:))
   disp('Sampled probability for the compressed qubits in the z-down direction is:')
   disp(probz2(1,:))
   disp(' ')
   disp(' ')
   disp('Sampled probability for the uncompressed qubits in the z-up direction is:')
   disp(probz1(2,:))
   disp('Sampled probability for the compressed qubits in the z-up direction is:')
   disp(probz2(2,:))
   disp(' ')
   disp(' ')
   disp('Sampled probability for the uncompressed qubits in the x-down direction is:')
   disp(probx1(1,:))
   disp('Sampled probability for the compressed qubits in the x-down direction is:')
   disp(probx2(1,:))
   disp(' ')
   disp(' ')
   disp('Sampled probability for the uncompressed qubits in the x-up direction is:')
   disp(probx1(2,:))
   disp('Sampled probability for the compressed qubits in the x-up direction is:')
   disp(probx2(2,:))
   disp(' ')
   disp(' ')
   disp('Sampled probability for the uncompressed qubits in the y-down direction is:')
   disp(proby1(1,:))
   disp('Sampled probability for the compressed qubits in the y-down direction is:')
   disp(proby2(1,:))
   disp(' ')
   disp(' ')
   disp('Sampled probability for the uncompressed qubits in the y-up direction is:')
   disp(proby1(2,:))
   disp('Sampled probability for the compressed qubits in the y-up direction is:')
   disp(proby2(2,:))
   disp(' ')
   disp(' ')
   disp('Number of samples for the compressed qubits in the z-direction is:')
   disp(shotcountz)
   disp('Number of samples for the compressed qubits in the x-direction is:')
   disp(shotcountx)
   disp('Number of samples for the compressed qubits in the y-direction is:')
   disp(shotcounty)
   disp('Number of samples for the uncompressed qubits is:')
   disp(shots1)
