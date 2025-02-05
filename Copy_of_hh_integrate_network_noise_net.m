%Written by Maxwell Sherman 6/2012
%Revised by Hyeyoung Shin 3/2014
%Revised by Stephanie Jones 3/2019

%Integrate the Hodgkin Huxley equations for a n cell network of excitory and
%inhibitory cells modeled with Il, Ina, Ik currents and GABAa/AMPA dynamic
%synapses. Includes Poisson noise in the I ce
function [tout,yout] = hh_integrate_network_noise_net(t, hh_param, syn_param, stim_param_i, stim_param_e,...
    connectivity_grid,pEI,pIE,pEE,pII,n,num_exc,num_inh,rate,noise_pop) % real connectivity is the connectivity (0,1) multiplied by probability

%%generate poisson spike times 
if noise_pop=='I'
    p_grid = poisson_spike_times_grid_icell(rate,300,num_inh); %generates array of poisson times for each i cell 
elseif noise_pop=='E'
    p_grid = poisson_spike_times_grid_ecell(rate,300,num_exc); %generates array of poisson times for each e cell
elseif noise_pop=='A'
    p_grid_a = poisson_spike_times_grid_ecell(rate,300,num_exc);
elseif noise_pop=='B'
    p_grid_b = poisson_spike_times_grid_ecell(rate,300,num_exc);
end

%following for loop steps through entries in spike time grid and adds each entry
%to a single vector, p_vector
p_vector_a = []; 
for i=1:length(p_grid_a)
    p_vector_a = [p_vector_a p_grid{i}];
end
p_vector_a = sort(p_vector_a); %sorts entries of p_vector in ascending order

p_vector_b = []; 
for i=1:length(p_grid_b)
    p_vector_b = [p_vector_b p_grid{i}];
end
p_vector_b = sort(p_vector_b); %sorts entries of p_vector in ascending order


% approx initial steady state for cells as initial conditions
X0 = zeros(1,6*n);
X0(1:n) = -80+rand(1,n)*20; %random initial condition for voltages
X0(n+1:n*2) = .0045+rand(1,n)*(.0403-.0045); %random initial conditions for n
X0(2*n+1:n*3) = .0010+rand(1,n)*(.0160-.0010); %random initial conditions for m
X0(3*n+1:n*4) = .9955+rand(1,n)*(.9970-.9955); %random initial conditions for h

%Load all parameter vectors outside of the integration function. Otherwise
%the vectors must be recreated each time the integration function iterates
%through (something on the order of hunreds of thousands of times). The
%parameter vectors become global variables and can are available for use
%within the nested integrator function. 

%parameters for E cells
GNa = hh_param(1);  
GK = hh_param(2);    
Gl = hh_param(3);   
vNa= hh_param(4); 
vK = hh_param(5);
vl = hh_param(6);

if pEI == 0
    gei=0;
else
    gei=syn_param(1)/(pEI*num_exc);
end
if pIE == 0
    gie=0;
else
    gie=syn_param(2)/(pIE*num_inh);
end
if pEE == 0
    gee=0;
else
    gee= syn_param(3)/(pEE*num_exc);
end
if pII == 0
    gii=0;
else
    gii = syn_param(4)/(pII*num_inh);
end
vampa=syn_param(5);
vgaba=syn_param(6);

%  applied currents, inhibitory and excitatory cells
I_I=stim_param_i(3);
I_E=stim_param_e(3);

%--------------------------------------------------------------------------

% Initialize variables used to calculate HH equations

%voltage variable indices (aka the index of the cells in the network stacked as a vector):
cell_index = zeros(n,1);
cell_index(1:n) = 1:n;
%n,m,h dynamic variable indices
n_index = zeros(n,1);
n_index(1:n) = n+1:n*2;
m_index = zeros(n,1);
m_index(1:n) = 2*n+1:n*3; 
h_index = zeros(n,1);
h_index(1:n) = 3*n+1:n*4;
%noise indices 
noise_index = zeros(n,1);
noise_index(1:n) = 5*n+1:n*6;

%Create vector that contains the synaptic terms for all cells--------------

%create gjk_grid that stores g value of each synapse
% n*n grid that's bidirectional 
gjk_grid=zeros(n,n);

gjk_grid(1:num_exc,1:num_exc)=gee;
gjk_grid(1:num_exc,num_exc+1:n)=gei;
gjk_grid(num_exc+1:n,1:num_exc)=gie;
gjk_grid(num_exc+1:n,num_exc+1:n)=gii;

    
%vsyn_vector has vgaba at indices for excitatory cells and
%vampa at indices for inhibitory cells.
vsyn_vector = zeros(n,1);
for i=1:num_exc
    vsyn_vector(i) = vampa; % set reversal potential of synapses
end
for i=num_exc+1:n
    vsyn_vector(i) = vgaba;
end

%convert to an n by n matrix (note every column is the same vector.
%the command repmat would accomplish the same goal, but is actually
%slower in this case. 
for i=1:n
    vsyn_grid(:,i)=vsyn_vector;
end
    
%create a n by n matrix where the ith column is a vector of entrie i 
current_cell=zeros(n,n);
for i=1:n
    current_cell(:,i)=ones(n,1)*i;
end

%create matrix of applied currents 
Iapp_matrix = zeros(n,1);
Iapp_matrix(1:num_exc) = I_E; %applied current to E cells
Iapp_matrix(num_exc+1:n) = I_I; % applied current to I cells

%for later use to load synaptic equations
e_cell_index = zeros(num_exc,1);
e_cell_index(1:num_exc)=1:num_exc;
    
i_cell_index = zeros(num_inh,1);
i_cell_index(1:num_inh)=num_exc+1:n; 


%=== begin nested functions ====

% Set up the Hodgkin-Huxley equations 
% Dynamic variables are stored in X = [v n m h]' %?

%Integration---------------------------------------------------------

%Use a for loop to break the tstep into time periods that correspond to the
%poisson spike times. We integrate from zero to the time of the first
% poisson spike. The integration is then stopped and restarted with
%the initial condition y0(i) = -pi, where is the index of the cell that has
%a poisson spike at that time(this forces a spike at the
%time of the poisson spike). The initial conditions for all other variables
%are set to their values at the end of the initial integration. Total time 
%and total y values are compounded in two vectors, tout and yout

tstart=0;
tout = []; %tout initially set to zero b/c time starts at 0
yout = []; %starts at whatever y0 is. In this case, 0

% I have 2 p vectors. How to integrate with them both? Loop thru current
% step time

for j = 1:length(p_vector)+1 % For each spike j, step ends at length(p)+1 so that integration 
                             % continues after the final spike
    if j<length(p_vector)+1 % If the integration is between poisson spikes
        [T Y] = ode23(@hhp,[tstart p_vector(j)],X0); % (1) solve hhp from tstart to the 
                                                     % current spike with IC X0
        % T is all the time for which there is a solution; Y is a matrix: [time, features]
        % where features are membrane potential - gating variables -
        % synaptic terms 
        %we not compile the produced data
        nt = length(T);
        tout = [tout;T(1:nt-1)]; 
        yout = [yout;Y(1:nt-1,:)]; % (2) append solution at all times except 
                                 % for the very last entry; becomes the 
                                 % first entry of the next iteration that
                                 % integrates with the next spike
    
        % (3) Reset IC X0
        % 1. First determines which cell is being forced to spike by searching
        % the vectors of the poisson array for the corresponding
        % spike time. 
        % 2. Then set corresponding noise value to 0.12 (maximal
        % conductance for ampa current)
        
       ind=[];
       for k=1:length(p_grid) %step through inhibitory cells %is k not time of spike???
            if any(p_vector(j)==p_grid{k})
                if noise_pop=='I' 
                 ind=[ind i_cell_index(k)]; % Indexing the i_cell_index to select neuron
                                            % indices that correspond to
                                            % time of spike 
                else
                 ind=[ind e_cell_index(k)]; 
                end
            end
        end
        
        X0=Y(nt,:);
        X0(5*n+ind)=.12; % why set this to 0.12 ?
        tstart = p_vector(j); % set tstart to j for the next iteration

    else %At the last spike index, integrate from the last spike to tstop
        [T Y] = ode23(@hhp,[tstart t],X0);
       
        nt = length(T);
        tout = [tout;T(1:nt)];
        yout = [yout;Y(1:nt,:)];
    end
end

function dXdt = hhp(t,X)
    dXdt = zeros(n*6,1); %preallocate dXdt vector
  
    %create a n by n matrix where every column are the Xs that
    %correspond to the synaptic gating variables.
    syn_terms = X(4*n+1:n*5);
    syn_grid=zeros(n,n);
    for l=1:n
        syn_grid(:,l)=syn_terms;
    end    
    
    %create the row vector containing synaptic terms for all cells within a
    %network 
    syn=sum(connectivity_grid.*gjk_grid.*syn_grid.*(X(current_cell)-vsyn_grid)...
        );
    % both between and within vsyn_grid are the same. They are of size
    % (n,n) for n neurons in pop A and n in pop B
    

%Create a matrix that contains dXdt for every cell except for the synaptic
%terms--------------------------------------------------------------------

    dXdt(1:n) = -GK*(X(n_index).^4).*(X(cell_index)-vK) - GNa*(X(m_index).^3).*X(h_index).*(X(cell_index)-vNa)...
        - Gl*(X(cell_index)-vl)+Iapp_matrix - X(noise_index).*X(cell_index);

%time for the unavoidable for loop to load the synaptic terms into each of
%dXdt equation------------------------------------------------------------
    
    for l = 1:n
        dXdt(l) = dXdt(l)-syn(l);
    end
    
% load equations for gating variables--------------------------------------
    
    dXdt(n+1:n*2) = an(X(cell_index)).*(1-X(n_index)) - bn(X(cell_index)).*X(n_index);

    dXdt(2*n+1:n*3) = am(X(cell_index)).*(1-X(m_index)) - bm(X(cell_index)).*X(m_index);

    dXdt(3*n+1:n*4) = ah(X(cell_index)).*(1-X(h_index)) - bh(X(cell_index)).*X(h_index);
    
    dXdt(5*n+1:n*6) = -.5*X(noise_index); %exponetial decay at a rate of tau=1/0.5=2;

%load synpatic variable terms----------------------------------------------    
    %some matrix manipulation here avoids for loops in loading synaptic gating equations.
    %A list of the E cell indices is created and fed through the kampa
    %function, and a list of I cell indices is created and fed through the
    %kgaba function

    dXdt(4*n+e_cell_index)=kampa(X(e_cell_index)).*(1-X(4*n+e_cell_index))-X(4*n+e_cell_index)/2;
    dXdt(4*n+i_cell_index)=kgaba(X(i_cell_index)).*(1-X(4*n+i_cell_index))-X(4*n+i_cell_index)/10;

return % returns the rate of change of X, which contains the 6 variables for all (n) cells 
end %hhp

%============
% Begin nested functions

%------------
% gating dynamics from Borgers et al 2008
%------------

%gatng terms for e cells
function val = an(v)  % opening rate functions of potassium activation
    val = 0.032*(v+52)./(1-exp(-0.2*(v+52)));
end

function val = bn(v)  % closing rate functions of potassium activation
    val = 0.5*exp(-0.025*(v+57));
end

function val = am(v)  % opening rate functions of sodium activation
    val = .32*(v+54)./(1-exp(-.25*(v+54)));
end

function val = bm(v)  % closing rate functions of sodium activation
    %SIGN TYPO IN PINTO JONES 2003, SIGN IS REVERSED
    val = (0.28*(v+27))./(exp(0.2*(v+27))-1);
end

function val = ah(v)  % opening rate functions of sodium inactivation
    val = 0.128*exp((-1/18)*(v+50));
end

function val = bh(v)  % closing rate functions of sodium inactivation
    val = 4./(1+exp(-0.2*(v+27)));
end

% synaptic terms

function val = kampa(v)  % opening rate functions of sodium inactivation
    val = 5*(1+tanh(v/4));
end

function val = kgaba(v)  % opening rate functions of sodium inactivation
    val = 2*(1+tanh(v/4));
end

%=== end nested functions ===
end %hh_integrate