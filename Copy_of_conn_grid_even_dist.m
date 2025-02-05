%Written by Maxwell Sherman 6/2012
%Revised by Hyeyoung Shin 3/2014

function connectivity_grid = conn_grid_even_dist(n,num_exc,probE_I,probI_E,probE_E,probI_I)

%This function returns the connectiviy matrix (connectivity_grid) of a network of n cells.
%The matrix is an n by n matrix.
%The row numbers correspond to pre-synaptic cell indices. The column
%numbers correspond to the post-synyaptic cell indices.
%Cells with index 1 through n_exc is excitatory, and the rest (index n_exc+1 through n) is inhibitory.
%The connectivity matrix is filled with ones and zeros. A 1 at position (i,j) means that there is a
%synapse between the ith and jth cells in the network.
%Ex. if there is a 1 at (1,2), cell #1 synapses on cell #2.

%This function produces evenly distributed connectivity. That is to say 
%each E cell will receive exactly the same number of I-E synapses 
%(probI_E*num_inh) and E-E synapses(probE_E*num_exc), and each I cell will 
%receive exactly probE_I*num_exc E_I synapses and probI_I*num_inh I_I synapses.

%-------------------------------------------------------------------------
%This function relies on a programing technique unique to MATLAB known as
%vectorizing. Rather than using for loops to repeat operations multiple
%times, a sinlge operation is performed an an entire vector or matrix all
%at the same time. This dramatically increases the speed and efficiency of
%the function.
%-------------------------------------------------------------------------

%Key:
    %n: total number of cells in network (excitatory + inhibitory)
    %num_exc: number of E cells in network
    %probE_I: probability of E to I connectivity
             %(i.e. (#of E_I synapses)/((#of E cells)*(#of I cells)) )
    %probI_E: probability of I to E connectivity
             %(i.e. (#of I_E synapses)/((#of I cells)*(#of E cells)) )
    %ProbE_E: probability of E to E connectivity
             %(i.e. (#of E_E synapses)/((#of E cells)*(#of E cells)) )
             %note, even though it is assumed that cells don't synapse on itself,
             %the denominator is chosen to be (#of E cells)*(#of E cells) for consistency.
             %same for I to I connections
    %ProbI_I: probability of I to I connectivity
             %(i.e. (#of I_I synapses)/((#of I cells)*(#of I cells)) )
    

%store cell layout as a vector. 1s are excitory, 0s are inhibitory
num_inh = n - num_exc;

%store connectivity as an n by n matrix 
connectivity_grid=zeros(n,n);

%create connectivity matrix
%connectivtiy for excitory neurons  
for i = 1:num_exc
    %i as the postsynaptic cell
    %I-E connections
    %for a given postsynaptic cell, randomly pick presynaptic inhibitory neuron
    ind = num_exc + randperm(num_inh, floor(probI_E*num_inh)); %in cell_vector, ihbitory cells have indices num_exc+1 <= index <= n 
                                      %round down #of synapses       
    connectivity_grid(ind,i)=1;
 
    %E_E connections
    %do not allow self connections
    ind = randi(num_exc, 1, floor(probE_E*num_exc)+1); %generate one extra random index in case i was included in list of presynaptic indices
    if find(ind==i) 
        ind=ind(ind~=i);
    else
        ind=ind(1:floor(probE_E*num_exc)); %if i wasn't included, take first floor(probE_E*num_exc) as index
    end
    connectivity_grid(ind,i) = 1;
end

%connectivtiy for inhibitory neurons  
for i = num_exc+1:n
    %i as the postsynaptic cell
    %E-I connections
    %for a given postsynaptic cell, randomly pick presynaptic excitatory neuron
    ind = randperm(num_exc, floor(probE_I*num_exc)); %round down #of synapses       
    connectivity_grid(ind,i)=1;
 
    %I_I connections
    %do not allow self connections
 
    ind = num_exc + randi(num_inh, 1, floor(probI_I*num_inh)+1); %generate one extra random index in case i was included in list of presynaptic indices
    if find(ind==i) 
        ind=ind(ind~=i);
    else
        ind=ind(1:floor(probI_I*num_inh)); %if i wasn't included, take first floor(probE_E*num_exc) as index
    end
    connectivity_grid(ind,i) = 1;
end

save connectivity_50cell_1_1_0_1.mat connectivity_grid n num_exc
