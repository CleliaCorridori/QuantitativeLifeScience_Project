%% Script for Phase diagram
clc
clear
close all

%patameters
sparsity=0.5;
epsilon=1e-5;
Npoints=100;
% Theoric values
dR2_vect = linspace(0,3, Npoints);
alpha_vect = linspace(1,7, Npoints);
psi = zeros(length(alpha_vect), length(dR2_vect));%initialization
for i=1:length(alpha_vect)
    for j=1:length(dR2_vect)
        dR2 = dR2_vect(j);
        par = theory(alpha_vect(i),epsilon,sparsity,dR2);
        psi(i,j)=par.psi;
    end
end

%% figure
figure(1);
logPsi = log(psi);
logPsi(logPsi<-10)=-Inf;
im = imagesc(dR2_vect([1,end]),alpha_vect([1,end]),logPsi);
%figure details
hold on
colorbar
axis xy
title('log \psi','FontSize',14)
ylabel('\alpha','FontSize',14)
xlabel('\delta R^2','FontSize',14)
text(0.5,5.5,'S','FontSize',18,'Color','k'); %S phase
text(2.5,3,'V','FontSize',18,'Color','k'); %V phase

% critical line
lambda=0:0.001:3;
alphaLine = 1./EE(lambda);
sparsity=0.5;
dR2Line = lambda.^2./(1-II(lambda)./EE(lambda))*(1-sparsity)/sparsity;
hold on;
plot(dR2Line,alphaLine,'w-','LineWidth',2.5)
