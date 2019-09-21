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

%Discretization for projection and some related one-time calculations
projsize = 20*xgridsize; %this doesn't have to be the same size as xgrid, 
%               but it's not meaningful to make it larger when using linear
%               interpolation/monotonicity constraint (since that's
%               determined by the values at the xgrid points; in between
%               the derivative sign is constant)
%               Also, changing is untested, may cause errors
projpoints = linspace(xgrid(1),xgrid(end),projsize);
proji = projpoints' * ones(1,projsize);
pdiff = abs(proji - proji');
%Tau regulates importance of closeness of function itself vs. derivative
%1 = equal weight
%For details see Delecroix et al. 1996
tau = .5;
H = ((tau^3)/2) * exp(-tau * pdiff / sqrt(2)) .* cos(tau * pdiff / sqrt(2) + pi/4);

%stepsize for derivative calculations (supports uneven grids)
diffx = diff(xgrid);

%define fitness function
fun = @(x) -abs(model2(gamma,epsilon,alpha,rho,sigmaerror,mu,IYssratio,KYssratio,nss,kss,delta,beta,b,xgrid,x) - .745);

%algorithm hyperparameters
ballsize = 1000000; %Maximum allowed perturbation.
                    %This value looks huge, but the units necessarily meaningful.
                    %Empirically, the perturbations this generates looks reasonable.
stepsize = .3;
mutationprob = .5;
popsize = 20;
generations = 10;

x = ugridbase + stepsize*randn(popsize,xgridsize);
%x = ugridbase;
tempxder = zeros(1,xgridsize);

%plot(x(1,:),'blue');
hold on
for p = 1:20
    %finite difference derivative approximation
    %"forward/backward on the ends, central in the middle"
    tempdiff =  diff(x(p,:));
    tempxder(1) = tempdiff(1) / diffx(1);
    tempxder(end) = tempdiff(end) / diffx(end);
    tempxder(2:end - 1) = (tempdiff(2:end) + tempdiff(1:end - 1)) ./ (diffx(1:end - 1) + diffx(2:end)); 
    xder = interp1(xgrid,tempxder,projpoints); %scale up the derivative to the right number of points
    
    if any(xder < 0) %project if it's decreasing somewhere
        alpha = quadprog(H,xder,[],[],[],[],zeros(projsize,1),[],[],options2);
        x(p,:) = x(p,:) + projdiff(tau,alpha,projpoints,xgrid);
    end
    %plot(x(p,:),'red');
    
    %recalculate derivative
    tempdiff =  diff(x(p,:));
    tempxder(1) = tempdiff(1) / diffx(1);
    tempxder(end) = tempdiff(end) / diffx(end);
    tempxder(2:end - 1) = (tempdiff(2:end) + tempdiff(1:end - 1)) ./ (diffx(1:end - 1) + diffx(2:end));
    norm = mean( (x(p,:) - ugridbase).^2 ) + (1/tau^2) * mean( (tempxder - 1./xgrid).^2 );
    if any(diff(x(p,:)) < 0)
        disp('monotonicity violated')
    end
    if norm > ballsize
        x(p,:) = ugridbase + ballsize * (x(p,:) - ugridbase)./norm;
    end
    plot(x(p,:))
end
        
% 
% %bernp = rand(popsize,1) <= mutationprob;
% %x = beta0*ones(popsize,1) + bernp .* (  rand(popsize,1)*( (ballhigh - balllow)*stepsize ) - balllow*stepsize  ) ;
% % %x = beta0*ones(popsize,1) + (  rand(popsize,1)*( (ballhigh - balllow)*stepsize ) - balllow*stepsize  ) ;
% 
% %x = linspace(balllow + .05,ballhigh - .05,popsize)' + (  rand(popsize,1)*( (ballhigh - balllow)*stepsize ) - balllow*stepsize  ) ;
% x = linspace(.8,.95,popsize)' + stepsize*(  (2*rand(popsize,1) - 1)*(ballhigh - balllow)  ) ;
% %x = linspace(balllow + .05,ballhigh - .05,popsize)' + (  randn(popsize,1) * sqrt(stepsize) ) ;
% F = zeros(popsize,generations);
% children = zeros(popsize - 1,1);
% parents = zeros(popsize - 1,1);
% xsave = zeros(popsize,generations); %for analyzing output
% for g = 1:generations
%     x(x < balllow) = balllow;
%     x(x > ballhigh) = ballhigh;
%     xsave(:,g) = x;   
%     parfor p = 1:popsize
%         F(p,g) = -fun(x(p));
%     end
%     [Fmax,Fargmax] = max(F(:,g));
%     %Fmax
%     g
%     %x(Fargmax)
%     probs = F(:,g) ./ sum(F(:,g));
%     parents = randsrc(popsize - 1,2,[linspace(1,popsize,popsize) ; probs']);
%     children(1:popsize - 1) = ( x(parents(:,1)).*F(parents(:,1),g) + x(parents(:,2)).*F(parents(:,2),g) ) ./ ( F(parents(:,1),g) + F(parents(:,2),g) );
%     mutate = rand(popsize - 1,1) <= mutationprob;    
%     children = children + mutate .* stepsize .*(  (2*rand(popsize - 1,1) - 1)*(ballhigh - balllow)  ) ;
%     %children = children + mutate .* (  randn(popsize - 1,1) * sqrt(stepsize) ) ;
%     x(1) = x(Fargmax);
%     x(2:popsize) = children;
% end
% 
% elitefitness = max(F,[],1);
% save('pop_and_fitness.mat', 'xsave','elitefitness')
 
toc