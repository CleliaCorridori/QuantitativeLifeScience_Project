% Simulation to compute the survivors species

clc
clear
close all

% Setting parameters
N=100;
dR2=1;
% dR2=3; It is possible to change the variance of the total supply
p=0.5; %sparsity

%% Decomment the desired total suppy (capacity)
%%% 1. R as in paper
capacity = [ones(N/2,1); -ones(N/2,1)].*sqrt(dR2/N)+1; %R

%%% 2. gaussian shape R
% x=linspace(1,50,50);
% base=linspace(-5,5,N);
% pd = makedist('Normal');
% norm_pdf = pdf(pd,base);
% capacity=0.9+norm_pdf;

%%% 3. sampled from gaussian distribution
% capacity = normrnd(1,1/N,[1,N]);

% dR2=std(capacity)^2;

%% Decomment the desired epsilon
epsilon_list = [1e-4, 0.03, 0.1]; %as in the paper
% epsilon_list = [0]; %to inspect P>N

%% theory results
alpha_list_theo = linspace(1,7,100);
survivors_theo = zeros(length(epsilon_list),length(alpha_list_theo));
for i=1:length(epsilon_list)
    epsilon = epsilon_list(i);
    for j=1:length(alpha_list_theo)
        alpha = alpha_list_theo(j);
        parameters = theory(alpha,epsilon,p,dR2);
        survivors_theo(i,j) = parameters.survivor_frac;
    end
end
data.survivors_theo = survivors_theo;
data.alphas_theo = alpha_list_theo;

%% simulation results
trials = 500; %500
alpha_list_sim = linspace(2,6,20);
%initialization of availability of resources (at equilibrium) and fraction
%of survivors:
availability_matx = zeros(length(epsilon_list),length(alpha_list_sim), trials, N); %initialization of h_i
survivors_sim = zeros(length(epsilon_list),length(alpha_list_sim),trials);

for i=1:length(epsilon_list)
    for j=1:length(alpha_list_sim)
        for s=1:trials
            epsilon=epsilon_list(i);
            alpha=alpha_list_sim(j);
            fprintf('epsilon = %f, alpha = %f, seed %d of %d.\n',epsilon_list(i), alpha_list_sim(j),s, trials);
            seed= round(1/(epsilon+alpha)+s); %recomputing the seed for the simulation
            org = Organisms_gen(N,alpha,epsilon,p,seed);
            org.capacity=capacity; %R
            [availability,~] = Equilibrium(org);
            availability_matx(i, j, s,:) = availability;
            delta = (org.enzMetStrat*availability(:) - org.cost);
            
            % Record the number of survivors
            delta_sorted = sort(delta, 'descend');
            surv = delta >= -1e-8; %Not zero, considering a certain threshold
            survivors_sim(i,j,s) = sum(surv); %per epsilon, alpha and trial fixed
        end
    end
end
data.survSim = survivors_sim/N; %fraction of survivors
data.alphas_sim = alpha_list_sim;
data.epsilon = epsilon_list;

%% figure
figure(5)
colors='brkg';
hold on

for i=1:length(epsilon_list) %cycling on epsilon
    %theory plot
    plot(data.alphas_theo,data.survivors_theo (i,:),'-','color',colors(i))
    
    % simulation plot
    mean_survivors = squeeze(mean(data.survSim(i,:,:),3));
    std_survivors = squeeze(std(data.survSim(i,:,:),[],3));
    ste = std_survivors./sqrt(size(data.survSim,3));
    errorbar(alpha_list_sim,mean_survivors,ste,'-','color',colors(i))
    
end

% other things for the plot:
% legend('1 theo','1 sim', '2 theo','2 sim','3 theo','3 sim','resources');
title('Fraction of survivors at equilibrium')
xlabel('\alpha')
ylabel('Fraction of survivors (P/N)')
yline(1,'-k')
grid on
% axis([2 6 0.55 1.05])
text(5.5,0.82,'\epsilon=0.1')
text(5.5,0.9,'\epsilon=0.03')
text(5.5,0.95,'\epsilon=10^{-4}')
% xline(3.7) %critical alpha

% For the epsilon proposed in the paper
% text(5.5,0.86,'\epsilon=0.1')
% text(5.5,0.94,'\epsilon=0.03')
% text(5.5,0.99,'\epsilon=10^{-4}')

% For epsilon=0
x=linspace(2,6,30);
hold on 
plot(x,x,'k-') %bisetrice
text(5.5,5,'\epsilon=0')
legend('data','y=x')
