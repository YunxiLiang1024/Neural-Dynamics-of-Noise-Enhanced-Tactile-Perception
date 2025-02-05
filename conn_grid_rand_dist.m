%Written by Maxwell Sherman 6/2012
%Revised by Hyeyoung Shin 3/2014

function connectivity_grid = conn_grid_rand_dist(n_each,num_exc,probE_I_within,probI_E_within,probE_E_across,probI_I_across)

%This function returns the connectiviy matrix (connectivity_grid) of a network of n cells.
%The matrix is an n by n matrix.
%The row numbers correspond to pre-synaptic cell indices. The column
%numbers correspond to the post-synyaptic cell indices.
%Cells with index 1 through n_exc is excitatory, and the rest (index n_exc+1 through n) is inhibitory.
%The connectivity matrix is filled with ones and zeros. A 1 at position (i,j) means that there is a
%synapse between the ith and jth cells in the network.
%Ex. if there is a 1 at (1,2), cell #1 synapses on cell #2.

%This function produces randomly distributed connectivity. That is to say,
%each E cell will not necessarily receive the same number of connections as
%every other E cell. Rather, the average probability of I-E connectivity
%across all E cells will be probI_E*num_inh while the average probability
%for E-E connectivity across all E cells will be probE_E*num_exc. The same
%pattern holds true for I cells. 

%-------------------------------------------------------------------------
%This function relies on a programing technique unique to MATLAB known as
%vectorizing. Rather than using for loops to repeat operations multiple
%times, a sinlge operation is performed an an entire vector or matrix all
%at the same time. This dramatically increases the speed and efficiency of
%the function.
%-------------------------------------------------------------------------

%Key:
    %r: number of rows in network
    %c: number of columns in network
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
             

num_inh = n_each - num_exc;

%store connectivity as an n by n matrix 
%preallocate matrix with 0s
connectivity_grid=zeros(2*n_each,2*n_each);

%set connection counts to zero
E_I_count=0;
I_E_count=0;
E_E_count=0;
I_I_count=0;

while E_I_count<floor(num_exc*probE_I_within)*num_inh
    %floor(num_exc*probE_I)*num_inh rather than floor(num_exc*probE_I*num_inh) to be consistent with conn_grid_even_dist
    %if floor(num_exc*probE_I*num_inh) != 0 (mod num_inh),
    %floor(num_exc*probE_I*num_inh) will end up imposing more synapses than conn_grid_even_dist.m
    %to prevent this, use floor(num_exc*probE_I)*num_inh instead
    %same applies for other types (I_E, E_E and I_I) of synapses
    a=randi(num_exc); %random presynaptic excitatory cell index
    b=randi([num_exc+1 n]); %random postsynaptic inhibitory cell index
    if connectivity_grid(a,b)==0 %check that this combination wasn't already connected
        connectivity_grid(a,b)=1;
        E_I_count = E_I_count+1;
    end
end

while I_E_count<floor(num_inh*probI_E)*num_exc
    a=randi([num_exc+1 n]); %random presynaptic inhibitory cell index
    b=randi(num_exc); %random postsynaptic excitatory cell index
    if connectivity_grid(a,b)==0 %check that this combination wasn't already connected
        connectivity_grid(a,b)=1;
        I_E_count = I_E_count+1;
    end
end

while E_E_count<floor(num_exc*probE_E)*num_exc
    a=randi(num_exc); %random presynaptic excitatory cell index
    b=randi(num_exc); %random postsynaptic excitatory cell index
    if a~=b && connectivity_grid(a,b)==0 %don't allow self connection & check that this combination wasn't already connected
        connectivity_grid(a,b)=1;
        E_E_count = E_E_count+1;
    end
end

while I_I_count<floor(num_inh*probI_I)*num_inh
    a=randi([num_exc+1 n]); %random presynaptic inhibitory cell index
    b=randi([num_exc+1 n]); %random postsynaptic inhibitory cell index
    if a~=b && connectivity_grid(a,b)==0 %don't allow self connection & check that this combination wasn't already connected
        connectivity_grid(a,b)=1;
        I_I_count = I_I_count+1;   
    end
end

% save connectivity_50cellR_1_1_0_0.mat connectivity_grid n num_exc
save connectivity_100cellR_1_1_0_0.mat connectivity_grid n num_exc


