function [Harvest, Delta] = Equilibrium(org)
% Optimization to compute equilibrium h_i using the Matlab function fmincon 

R = (org.capacity(:))'; %resources
N = org.N; %number of resources
fun = @(h)funToMinimize(h, R);
hessianfcn = @(h,lambda)Hessian(h, R);

tolVal = 1e-10;
%set the option to optimize the fmincon solver
%help at: https://www.mathworks.com/help/optim/ug/optimization-options-reference.html
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'HessFcn', hessianfcn, ...
    'ConstraintTolerance',tolVal,'OptimalityTolerance',tolVal,'StepTolerance',tolVal,'Display','off');
h0 = 0.9*ones(1,N);
lb = zeros(1,N);
ub = Inf(1,N);
step_rank = 5000;

met_strat = double(org.enzMetStrat);
costs = org.cost;
% sort by increasing cost
[costs, order] = sort(costs,'ascend');
met_strat = met_strat(order,:);
maxRank = min(step_rank, length(order)); %select the number of species considered
done = false;
while ~done
    subset = 1:maxRank;
    % when optimizing, use only a subset of species:
    %"fmincon": Find minimum of constrained nonlinear multivariable function
    Harvest = fmincon(fun,h0,met_strat(subset,:),costs(subset),[],[],lb,ub,[],options);
    Delta = met_strat*Harvest'-costs; %surplus
    
    % verify that ALL the constraints are satisfied
    done = all(Delta<=0);
    if ~done
        fprintf('Some ignored constraints are violated; re-trying.\n');
        maxRank = min(maxRank+step_rank,length(order)); % increase the subset
    end
end 
Delta(order) = Delta;
end


function [f,grad] = funToMinimize(h, R)
% R Capacity
f = -R*log(h(:));

grad = -R(:)./h(:); %gradient
end

function H = Hessian(h,R)
    H = diag(R(:)./(h(:).^2));
end
