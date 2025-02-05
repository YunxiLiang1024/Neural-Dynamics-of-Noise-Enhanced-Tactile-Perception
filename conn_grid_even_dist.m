%Written by Maxwell Sherman 6/2012
%Revised by Hyeyoung Shin 3/2014

function connectivity_grid = conn_grid_even_dist(n_each,num_exc,connectivity_scaling_within,connectivity_scaling_between)

rng(22);
pAEAI=connectivity_scaling_within(1);
pAIAE=connectivity_scaling_within(2);
pBEBI=connectivity_scaling_within(3);
pBIBE=connectivity_scaling_within(4);

pAEBI=connectivity_scaling_between(1);
pAEBE=connectivity_scaling_between(2);
pAIBE=connectivity_scaling_between(3);
pAIBI=connectivity_scaling_between(4);
pBEAI=connectivity_scaling_between(5);
pBEAE=connectivity_scaling_between(6);
pBIAE=connectivity_scaling_between(7);
pBIAI=connectivity_scaling_between(8);


% implements a network where there is EI/IE within network, and purely EE
% II between networks. Disconnect the EE, II within a network 

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
    
% Make the first 50 cells network A and the last 50 cells network B

%store cell layout as a vector. 1s are excitory, 0s are inhibitory

%store connectivity as an n by n matrix 
% 1-40 Aexc 41-50 Ainh; 51-90 Aexc 91-100 Ainh
connectivity_grid=zeros(n_each*2,n_each*2);
num_inh=n_each-num_exc;

% i is postsynaptic, ind is presynaptic
% -------------------AA connections----------------
for i = 1:num_exc
    % AI-AE
    ind = num_exc + randperm(num_inh, floor(pAIAE*num_inh));       
    connectivity_grid(ind,i)=1; 
end
for i = num_exc+1:n_each %target is 41-50
    %AE-AI
    ind = randperm(num_exc, floor(pAEAI*num_exc)); %round down #of synapses       
    connectivity_grid(ind,i)=1;
end
% -------------------BB connections----------------
for i = n_each+1:n_each+num_exc %51-90
    % BI-BE 
    ind = n_each+num_exc + randperm(num_inh, floor(pBIBE*num_inh)); 
    connectivity_grid(ind,i)=1; % 1 indicates there's a connection  
end
for i = n_each+num_exc+1:n_each*2
    %E-I connections in pop B
    %for a given postsynaptic cell, randomly pick presynaptic excitatory neuron
    ind = n_each+randperm(num_exc, floor(pBEBI*num_exc)); %round down #of synapses       
    connectivity_grid(ind,i)=1;
end
% ---------------B->A connections-----------------
% BE - AE
for i = 1:num_exc
    ind = randi([n_each+1,n_each+num_exc], floor(pBEAE*num_exc)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
% BE - AI
for i = num_exc+1:n_each
    ind = randi([n_each+1,n_each+num_exc], floor(pBEAI*num_exc)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
% BI - AE
for i = 1:num_exc
    ind = randi([n_each+num_exc+1,2*n_each], floor(pBIAE*num_inh)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
% BI - AI
for i = num_exc+1:n_each
    ind = randi([n_each+num_exc+1,2*n_each], floor(pBIAI*num_inh)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
% ---------------A->B connections-----------------
% AE - BE
for i = n_each+1:n_each+num_exc
    ind = randi([1,num_exc], floor(pAEBE*num_exc)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
% AE - BI
for i = n_each+num_exc+1:2*n_each
    ind = randi([1,num_exc], floor(pAEBI*num_exc)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
% AI - BE
for i = n_each+1:n_each+num_exc
    ind = randi([num_exc+1,n_each], floor(pAIBE*num_inh)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
% AI - BI
for i = n_each+num_exc+1:2*n_each
    ind = randi([num_exc+1,n_each], floor(pAIBI*num_inh)); % select exc indices in popB
    connectivity_grid(ind,i) = 1;
end
save connectivity_100cell_1_1_mixed1.mat connectivity_grid n_each num_exc
