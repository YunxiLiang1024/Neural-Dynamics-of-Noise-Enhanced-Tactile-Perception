function noisy_spike_times_grid = poisson_spike_times_grid_icell(rate,time,num_inh)
%function generates a cell array of spike times based on a random poisson
%distribution for inibitory neurons in the system.
%modified by Stephanie Jones 3/2019

%rate is the expected rate of noise-driven spikes in our system
%(note: rate should be reported in hz)
%Time is the duration of our model in ms.



%convert rate into ISI (ISI is the average interspike time thatwill 
%be used to generate a vector of interspike times  

lamda=rate; %lamda is number of cylces per second
ISI = 1/lamda; %computes average interspike time in seconds
ISI = ISI*1000; %converts IST to ms


noisy_spike_times_grid = cell(num_inh,1); %preallocate cell array to store vectors of spike times

%use a for loop to step through cell grid. If entry is 1 (i.e. if position
%is an excitatory cell), generate vector of spike times. Else, cell must be
%inhibitory and not subject to noise

%for loop through inhibitory cells.
%a while loop is used to sum interspike times generated from a
%poisson distribution. Interspike times are generated while the their sum is less than the time window.
%the running sum is stored as an entry in the vector noisy_spike_times
%where the index is equivalent to the number of times the loop has been
%run through. This vector is then stored in the array noisy_spike_times_gird

for j = 1:num_inh
    count = exprnd(ISI); %first spike time
    %exprnd generates possion distributed interspike times
    i=1; %initialize indexing variable
    noisy_spike_times = [];
    while count<time
        noisy_spike_times(i) = count;
        i=i+1;
        count = count+exprnd(ISI); %exprnd generates possion distributed interspike times
    end
    noisy_spike_times_grid{j} = noisy_spike_times;
end
