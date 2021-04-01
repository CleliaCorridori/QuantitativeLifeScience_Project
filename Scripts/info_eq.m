function [availability, deltas] = info_eq(N,alpha, epsilon, capacity, trials_tot)
%Function to compute the resource availability and the delta of survivors,
%used in Simulation_script

% initialize:
availability = zeros(trials_tot, N);
deltas=zeros(trials_tot,round(N*alpha));
for trial=1:trials_tot
    fprintf('alpha = %f, epsilon = %f.4, seed %d of %d.\n',alpha, epsilon, trial, trials_tot);
    %computation h_i at equilibrium
    seed= round(1/(epsilon+alpha)+trial); %as done before for the previous figure, change the initial seed
    org = Organisms_gen(N,alpha,epsilon,0.5,seed); %0.5 sparsity
    org.capacity=capacity; %R
    h = Equilibrium(org);
    availability(trial,:) = h; %equilibrium point, for each trial
    
    % Resource Surplus:
    delta = (org.enzMetStrat*h(:) - org.cost); %%%
    [delta_sorted,~] = sort(delta, 'descend'); %cost for each species, sorted
    deltas(trial,:) = delta_sorted;
end

end