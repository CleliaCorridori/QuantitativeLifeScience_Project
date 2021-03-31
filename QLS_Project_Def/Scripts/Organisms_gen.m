% Simulation of an Ecosystem with many species
function org = Organisms_gen(N,alpha,epsilon,sparsity,seed)
    rng(seed);
    org.N = N;
    org.P = round(alpha*N);
    org.enzMetStrat = met_strat_gen(N,org.P,sparsity); %metabolic strategies Sigma:_mu
    
    org.enzTotSpec = sum(org.enzMetStrat,2); %number of used resource per species
    
    for i=1:length(org.enzTotSpec) %met strat normalization
        org.enzMetStrat_norm(i,:)=org.enzMetStrat(i,:)/org.enzTotSpec(i);
    end
    org.enzMetStrat=org.enzMetStrat_norm;
    
    org.cost = Costs_gen(org, epsilon);  %Chi, Cost
    org.cost=org.cost./org.enzTotSpec;
end

%% metabolic strategies
function met_strat_matx = met_strat_gen(N,P,sparsity) 
    % Organisms: each pathway with probability 1/2
    met_strat_matx = rand(P,N)<sparsity;
    
    % Check zeros: if Yes, recompute the function
    zero = any(sum(met_strat_matx,2)==0);
    if zero
        fprintf('\t Empty organism generated. Retrying.\n');
        met_strat_matx = met_strat_gen(params);
    end
end

%% chi: cost
function costs = Costs_gen(org, epsilon)
    x_rand = random('Normal',0,1,[org.P,1]);
    costs = org.enzTotSpec - epsilon*x_rand;
    if any(costs<0)
        fprintf('Negative costs generated. Repeating.\n');
        costs = Costs_gen(org, epsilon);
    end
end


