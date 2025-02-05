function noisy_spike_times_grid = poisson_spike_times_grid_ecell(rate,time,n_each)

%Written by Maxwell Sherman 6/2012
%Revised by Hyeyoung Shin 3/14
%Revised by Stehanie Jones 3/2019

%function generates a cell array of spike times based on a random poisson
%distribution for excitory neurons in the system.

%rate is the expected rate of noise-driven spikes in our system
%(note: rate should be reported in hz)
%Time is the duration of our model in ms.

%--------------------------------------------------------------------------
%What is a Poisson distribution: a Poisson distribution is a distribution
%that expresses the probability of a given number of events ocurring in a
%fixed interval or time and/or space. The distribution requires three
%assumptions. 1) Events occur at a known average rate. 2)Events occur
%independently of each other. 3) The probability of an event occuring
%becomes exponentially less likely as one moves away from the average. 

%How does it relate to neural modeling: a simple way of modeling a spike
%train for a cell is to specify the expected average rate of firing and
%then add stochastic variation to it via a Poisson distribution. This
%means that the probability that the cell fires more or less than the
%average decreases exponentially as the absolute value of the difference
%between the average and event number increases linearly. 
%--------------------------------------------------------------------------


%convert rate into lambda (lambda is the average interspike time that will 
%be used to generate a vector of interspike times)  

lamda=rate; %lamda is number of cylces per second
ISI = 1/lamda; %computes average interspike time in seconds
ISI = ISI*1000; %converts IST to ms


noisy_spike_times_grid = cell(n_each,1); %preallocate cell array to store vectors of spike times
%each row of cell array will store a vector (variable length for each excitatory cell)
%of noisy_spike_times for that excitatory cell

%for loop through excitatory cells.
%a while loop is used to sum interspike times generated from a
%poisson distribution. Interspike times are generated while the their sum is less than the time window.
%the running sum is stored as an entry in the vector noisy_spike_times
%where the index is equivalent to the number of times the loop has been
%run through. This vector is then stored in the array noisy_spike_times_gird
for j = 1:n_each
    count = exprnd(ISI); %first spike time, sampled randomly from an exponential distribution with mean ISI
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
