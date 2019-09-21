clear
close all
rng(12378)
%parpool(20)

tic
options = optimset('Display','iter');
options2 = optimset('Display', 'off');

%This set of parameters (gamma through b) from Fabrizio Zilibotti homework
gamma       = 1;                %assuming log utility base. Modification currently unsupported.
epsilon     = 2;                %elasticity of labor supply
alpha       = .4;               %capital's share of income
rho         = .95;
sigmaerror  =.01;
mu          = 0;                %long-run growth rate
%Steady State Values
IYssratio   = .25;
KYssratio   = 14;
nss         = 1/3;
kss         = 14^(.06) / 3;              
%Calibrated Parameter Values
delta   = 1/56;                  %depreciation rate
beta    = .964;                  %discount rate
b       = 7.2;                  %disutility of labor multiple

%Discretization for the utility function
xgridsize = 50;
xgrid = linspace(.00005,1.0,xgridsize);
ugridbase = log(xgrid); %initial value = log utility

%define fitness function
fun = @(x) -abs(model2(gamma,epsilon,alpha,rho,sigmaerror,mu,IYssratio,KYssratio,nss,kss,delta,beta,b,xgrid,x) - .745);

%algorithm hyperparameters
ballsize = 1; %Maximum allowed overall perturbation.
stepsize = .2;
popsize = 20;

perturbation = zeros(popsize,xgridsize - 1);
increase = diff(ugridbase);

for p = 1:popsize
    for index = 1:xgridsize - 1
        candidate = stepsize*(2*rand(1) - 1);
        maxdecrease = increase(index);
        if stepsize < maxdecrease
            perturbation(p,index) = stepsize*candidate;
        else
            perturbation(p,index) = maxdecrease*candidate;
        end        
    end
    norm = sqrt(mean( (perturbation(p,:) ./ increase ).^2));
    if norm > ballsize       
        perturbation(p,:) = ballsize * perturbation(p,:)/norm;
    end
end

x = zeros(popsize,xgridsize);
x(:,1) = ugridbase(1);
for index = 2:xgridsize
    x(:,index) = x(:,index - 1) + increase(index - 1) + perturbation(:,index - 1);
end

plot(perturbation')
figure(2)
plot(perturbation(:,5:end)')
figure(3)
plot(x')
figure(4)
plot(x(:,5:end)')
 
toc