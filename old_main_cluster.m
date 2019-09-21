rng(12378)
parpool('cluster1',20)

tic
options = optimset('Display','iter');

gamma       = 1;                %assuming log utility
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

fun = @(x) -abs(model(gamma,epsilon,alpha,rho,sigmaerror,mu,IYssratio,KYssratio,nss,kss,delta,x,b) - .745);
%[test, argtest] = fminbnd(fun,.8,.999,options)
%test = model(gamma,epsilon,alpha,rho,sigmaerror,mu,IYssratio,KYssratio,nss,kss,delta,beta,b)

balllow = .5;
ballhigh = .99;
stepsize = .1;
mutationprob = .5;
popsize = 20;
beta0 = .964;
generations = 10;

%bernp = rand(popsize,1) <= mutationprob;
%x = beta0*ones(popsize,1) + bernp .* (  rand(popsize,1)*( (ballhigh - balllow)*stepsize ) - balllow*stepsize  ) ;
% %x = beta0*ones(popsize,1) + (  rand(popsize,1)*( (ballhigh - balllow)*stepsize ) - balllow*stepsize  ) ;

%x = linspace(balllow + .05,ballhigh - .05,popsize)' + (  rand(popsize,1)*( (ballhigh - balllow)*stepsize ) - balllow*stepsize  ) ;
x = linspace(.8,.95,popsize)' + stepsize*(  (2*rand(popsize,1) - 1)*(ballhigh - balllow)  ) ;
%x = linspace(balllow + .05,ballhigh - .05,popsize)' + (  randn(popsize,1) * sqrt(stepsize) ) ;
F = zeros(popsize,generations);
children = zeros(popsize - 1,1);
parents = zeros(popsize - 1,1);
xsave = zeros(popsize,generations); %for analyzing output
for g = 1:generations
    x(x < balllow) = balllow;
    x(x > ballhigh) = ballhigh;
    xsave(:,g) = x;   
    parfor p = 1:popsize
        F(p,g) = -fun(x(p));
    end
    [Fmax,Fargmax] = max(F(:,g));
    %Fmax
    g
    %x(Fargmax)
    probs = F(:,g) ./ sum(F(:,g));
    parents = randsrc(popsize - 1,2,[linspace(1,popsize,popsize) ; probs']);
    children(1:popsize - 1) = ( x(parents(:,1)).*F(parents(:,1),g) + x(parents(:,2)).*F(parents(:,2),g) ) ./ ( F(parents(:,1),g) + F(parents(:,2),g) );
    mutate = rand(popsize - 1,1) <= mutationprob;    
    children = children + mutate .* stepsize*(  (2*rand(popsize,1) - 1)*(ballhigh - balllow)  ) ;
    %children = children + mutate .* (  randn(popsize - 1,1) * sqrt(stepsize) ) ;
    x(1) = x(Fargmax);
    x(2:popsize) = children;
end

elitefitness = max(F,[],1);
save('pop_and_fitness.mat', 'xsave','elitefitness')

toc