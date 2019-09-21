clear
close all
rng(12378)
parpool(20)

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
ballsize = 1; %Maximum allowed perturbation.
stepsize = .2;
mutationprob = .2;
popsize = 20;
generations = 300;

perturbation = zeros(popsize,xgridsize - 1);
increase = diff(ugridbase);

%initial perturbations
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

%initialize the population
x = zeros(popsize,xgridsize);
x(:,1) = ugridbase(1);
for index = 2:xgridsize
    x(:,index) = x(:,index - 1) + increase(index - 1) + perturbation(:,index - 1);
end

%reduce the stepsize; we want the initial population spread out but don't
%want to move around too randomly during the rest of the algorithm
stepsize = .1;


F = zeros(popsize,generations);
children = zeros(popsize - 1,xgridsize);
parents = zeros(popsize - 1,xgridsize);
mutation = zeros(xgridsize - 1,1);
xsave = zeros(popsize,xgridsize,generations); %for analyzing output
for g = 1:generations
    xsave(:,:,g) = x;   
    parfor p = 1:popsize
        F(p,g) = -fun(x(p,:));
    end
    [Fmax,Fargmax] = max(F(:,g));
    %Fmax
    g
    %x(Fargmax)
    probs = F(:,g) ./ sum(F(:,g));
    parents = randsrc(popsize - 1,2,[linspace(1,popsize,popsize) ; probs']);
    children = ( x(parents(:,1),:).*F(parents(:,1),g) + x(parents(:,2),:).*F(parents(:,2),g) ) ./ ( F(parents(:,1),g) + F(parents(:,2),g) );    
    mutate = rand(popsize - 1,1) <= mutationprob;    
    if any(mutate > 0)        
        mutation_indices = find(mutate > 0); 
        for m = 1:size(mutation_indices)     
            pre_mutate_increase = diff(children(m,:));
            %generate mutation
            for index = 1:xgridsize - 1
                candidate = stepsize*(2*rand(1) - 1);
                maxdecrease = pre_mutate_increase(index);
                if stepsize < maxdecrease
                    mutation(index) = stepsize*candidate;
                else
                    mutation(index) = maxdecrease*candidate;
                end        
            end
            %apply it
            for index = 2:xgridsize
                children(m,index) = children(m,index - 1) + pre_mutate_increase(index - 1) + mutation(index - 1);
            end
        end
    end
            
    x(1,:) = x(Fargmax,:);
    x(2:popsize,:) = children;
    
    for n = 2:popsize
        perturbation_temp = increase - diff( x(n,:) );
        norm = sqrt(mean( (perturbation_temp ./ increase ).^2));
        if norm > ballsize
            %renormalize
            perturbation_temp(p,:) = ballsize * perturbation_temp(p,:)/norm;
            x(n,1) = ugridbase(1);
            for index = 2:xgridsize
                x(n,index) = x(:,index - 1) + increase(index - 1) + perturbation_temp(:,index - 1);
            end
        end
    end
    if mod(g,50) == 0  
        %checkpoint
        elitefitness = max(F,[],1);
        save('checkpoint_pop_and_fitness.mat', 'xsave','elitefitness')
    end
end

elitefitness = max(F,[],1);
save('pop_and_fitness.mat', 'xsave','elitefitness')
 
toc