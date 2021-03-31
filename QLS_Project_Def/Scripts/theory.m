function params = theory(alpha,epsilon,p,dR2) % p is the sparsity
%Function to compute of the theoretical values
x   = @(lambda) 1-alpha*EE(lambda); % eq S17 Supporting materials
psi = @(lambda) -2*(1-p)*(1-alpha*II(lambda))./(alpha*dII(lambda))-(1-p)*lambda; % eq S16 Supporting materials
eq_zero = @(lambda) sqrt(p*(1-p)./(1-alpha*II(lambda)).*(x(lambda).^2)*dR2 + epsilon^2./(1-alpha*II(lambda))) - psi(lambda); %equal to zero

%check on 1-alpha*II(lambda)>0
lambdaMin = 1e-10 + fzero(@(lambda)1-alpha*II(lambda),10);
assert(1-alpha*II(lambdaMin)>0,'Error!'); %this print an error, if it occurs

lambda0 = fzero(eq_zero,[lambdaMin,10]); %solution of equation eq_zero to find lambda

%%%%% output
params.lambda = lambda0;
params.x = x(lambda0);
params.I = II(lambda0);
params.psi = psi(lambda0);                      %psi
params.q = (params.psi^2-epsilon^2)/(p*(1-p));  %q
params.mStar = lambda0*params.psi/p;            %m
params.survivor_frac = alpha*EE(lambda0);       %fraction of surviviors

end

function x=EE(lambda)
% Function to compute the function signaled in the paper as E(lambda)
x = 1/2.*erfc(lambda/sqrt(2));
end
function x=II(lambda)
% Function to compute the function signaled in the paper as I(lambda)
x = - lambda.*exp(-0.5*lambda.^2)/sqrt(2*pi) + (1+lambda.^2)/2.*erfc(lambda/sqrt(2));
end
function x=dII(lambda)
% Function to compute the derivative I(lambda)
x = - exp(-0.5*lambda.^2)*sqrt(2/pi) + lambda.*erfc(lambda/sqrt(2));
end

