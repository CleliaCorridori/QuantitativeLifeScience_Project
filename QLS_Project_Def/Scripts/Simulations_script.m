% simulation figure 3:
clc
clear
close all

% Setting parameters
N=50;
% dR2=1;
% dR2=1; It is possible to change the variance of the total supply
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

%% Fixed total influx of resource
figure(1)
plot(capacity,'-*b')
grid on
xlabel('Resources')
ylabel('R_i')
title('External supply')
axis([0 N 0.7 1.3])  % set it
% yline(1,'k','LineWidth',2)
yline(mean(capacity),'k','LineWidth',1)

%% Phase transition or crossover for different epsilon, increasing alpha
% THEORY PART
val_alpha_teo = logspace(0,2,200); %alpha between 10^0 and 10^2, #points=200
alpha = 1; %initialization of alpha
max_exp=5; %exp of epsilon
psi=zeros(1,length(val_alpha_teo));
Psi=zeros(length(1:max_exp),length(val_alpha_teo));
for exp_eps = 1:max_exp 
    epsilon = 10^(-exp_eps); %loop on different values of epsilon
    for i=1:length(val_alpha_teo) %loop on different values of alpha
        alpha = val_alpha_teo(i);
        params(i) = theory(alpha,epsilon,p,dR2);
        psi(i)=params.psi;
    end
    Name_Par = sprintf('paramsList%d',exp_eps);
    data_diff_eps.(Name_Par) = params;
    Psi(exp_eps,:)=psi(:);
end
%% Figure: theory part
figure(2)
newcolors = [0, 0.4470, 0.7410
             0.8500, 0.3250, 0.0980
             0.9290, 0.6940, 0.1250
             0.4940, 0.1840, 0.5560
             0.4660, 0.6740, 0.1880];
%              0.3010, 0.7450, 0.9330]; %decomment for a new color (sesto)
colororder(newcolors)

loglog(val_alpha_teo,[data_diff_eps.paramsList1.psi],'-','LineWidth',2); %power 1
hold on
loglog(val_alpha_teo,[data_diff_eps.paramsList2.psi],'-','LineWidth',2); %power 2
loglog(val_alpha_teo,[data_diff_eps.paramsList3.psi],'-','LineWidth',2); %power 3
loglog(val_alpha_teo,[data_diff_eps.paramsList4.psi],'-','LineWidth',2); %power 4
loglog(val_alpha_teo,[data_diff_eps.paramsList5.psi],'-','LineWidth',2); %power 5

%text box
%-1
text(50,0.25,'\epsilon=0.1')
%-2
text(50,0.025,'\epsilon=0.01')
% -3
text(50,0.0025,'\epsilon=10^{-3}')
% -4
text(50,0.00025,'\epsilon=10^{-4}')
% -5
text(50,0.000025,'\epsilon=10^{-5}')

grid on
xlabel('\alpha')
ylabel('\Psi')
title('\Psi in function of \alpha for different \epsilon')

%% SIMULATION PART
val_alpha_sim = [1 2 3 4 6 8 10 20 50 100]; %values of alpha %5,15
N_alpha=length(val_alpha_sim);
repetitions = 10; % ### trials
availability_matx = zeros(max_exp,N_alpha, repetitions, N); %initialization of h_i
for alphaIdx = 1:N_alpha
    alpha=val_alpha_sim(alphaIdx);
    for trial = 1:repetitions
        fprintf('alpha = %d, trial %d of %d.\n',alpha,trial,repetitions);
        for exp_eps = 1:max_exp %loop on epsilon
            epsilon=10^(-exp_eps);
            seed= round((1/epsilon+alphaIdx)+trial); %setting the seed for the generation of Chi and sigma
            
            org = Organisms_gen(N,alpha,epsilon,p,seed);
            org.capacity=capacity; %R, total influx of resources
            
            availability = Equilibrium(org); % h
            availability_matx(exp_eps, alphaIdx, trial,:) = availability; %save h into a matrix
        end
    end
end

%% FIGURE for simulation part
for i=1:max_exp
    Name_Par = sprintf('paramsList%d',i);
    params = data_diff_eps.(Name_Par);

    alpha_sim = repmat(val_alpha_sim,[size(availability_matx,3),1])';
    m_i=1-availability_matx(i, :, :, :); %for each epsilon
    q_sim = squeeze(sum(m_i.^2,4)); %computing q, eq 3 of the paper
    epsilon = 10^(-i);
    psi_sim = sqrt(q_sim(:)*p*(1-p)+epsilon^2); %using relation with m,q and epsilon to compute Psi
    plot(alpha_sim(:), psi_sim(:),'.')
end

%% Equilibrium availability of resources h_i
% N=50
% mean and std over 500 trials
phaseVreplicas = 500;      % #repetiotions, equal for both phases
phaseSreplicas = 500;

selectAlpha = [2 10]; %for V and S phase, form paper
epsilon = 1e-3; %from paper 
availability_500.V= info_eq(N,selectAlpha(1), epsilon, capacity, phaseVreplicas); %V hi
[availability_500.S, deltas] = info_eq(N,selectAlpha(2), epsilon,capacity, phaseSreplicas); %S hi+Delta

%% Figure for equilibrium availability of resources h_i
figure(3)
Avail_V=availability_500.V; %V
Avail_S=availability_500.S; %S

% V phase data
%mean over trials
availability_mean_V = mean(Avail_V);
%index for different substrates of the resource
res_V = 1:length(availability_mean_V);
% std values
availability_std_V = std(Avail_V);
upper_V = availability_mean_V+availability_std_V;
lower_V = availability_mean_V-availability_std_V;
%plot
hold on
p1=plot(res_V,availability_mean_V,'b.','LineWidth',2);
patch([res_V, fliplr(res_V)],[upper_V, fliplr(lower_V)],'b','FaceAlpha',0.1,'EdgeColor','none');

%S phase data
%mean, std over trials
availability_mean_S = mean(Avail_S);
availability_std_S = std(Avail_S);
res_S = 1:length(availability_mean_S); %index for different substrates
upper_S = availability_mean_S+availability_std_S;
lower_S = availability_mean_S-availability_std_S;
%plot
patch([res_S, fliplr(res_S)],[upper_S, fliplr(lower_S)],'r','FaceAlpha',0.1,'EdgeColor','none');
p2=plot(res_S,availability_mean_S,'r.','LineWidth',2);
legend([p1 p2],{'V phase','S phase'})
grid on
xlabel('Resources')
ylabel('h_i')
title('Resource availability')

%% Distribution of the resource surplus: just extincted species
figure(4)
alpha_S = selectAlpha(2); % focus on S phase
epsilon = 1e-3; 

% theory
theory_data = theory(alpha_S, epsilon, p, dR2);
mean_theo = theory_data.lambda * theory_data.psi;
step = 0.04;
binmean = (-5:step:5)*mean_theo;
%normal distribution with:  mu = mean & sigma=psi
distrib = pdf('normal',binmean, -mean_theo, theory_data.psi)*step*mean_theo;

counts = histcounts(deltas(:),binmean)/(alpha_S*N*size(deltas,1));
plot(binmean, distrib,'k-','LineWidth',2)
hold on
plot(binmean(1:end-1)+diff(binmean),counts,'r*');
axis([-0.007,0,0,0.03]);
xline(-mean_theo,'b','LineWidth',2)
legend('theoretical distribution','data from simulation','mean','location','northwest')
title('p(\Delta)')
xlabel('\Delta')
ylabel('p(\Delta)')

