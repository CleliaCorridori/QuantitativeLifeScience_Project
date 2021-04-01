function x=II(lambda)
% Function to compute the function signaled in the paper as I(lambda)
x = - lambda.*exp(-0.5*lambda.^2)/sqrt(2*pi) + (1+lambda.^2)/2.*erfc(lambda/sqrt(2));
end