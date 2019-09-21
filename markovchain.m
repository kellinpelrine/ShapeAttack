function [state, transition, invdist] = markovchain(M, rho, sigma, suppwidth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% markovchain.m
% by Jen-Jen La'O
%
% creates a discrete M-State Markov chain to approximate
% a first-order autoregressive process:
% state_t = rho*state_(t-1) + u_t
% u_t is error term (white noise): iid N(0,sigma)
%
% Input to 'markovchain':
% 'M' is number of states of markov chain
% 'rho' is AR(1) coefficient
% 'sigma' is standard deviation of error term u_t
% 'suppwidth' is width of support of state space in terms of std devs
% Output of 'markovchain':
% 'state' is the 1xM state space vector
% 'transition' is the MxM transition matrix of the markov chain
% 'invdist' is the 1xM invariant distribution vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create State Space
stdev=sqrt((sigma^2)/(1-rho^2));
smax=suppwidth*stdev;
smin=-suppwidth*stdev;
interval=(smax-smin)/(M-1);

state=[smin:interval:smax];

%Calculate Transition Matrix
for j=1:M
    transition(j,1)=normcdf(state(1)-rho*state(j)+interval/2,0,sigma);
    for k=2:M-1
        transition(j,k)=normcdf(state(k)-rho*state(j)+interval/2,0,sigma)-normcdf(state(k)-rho*state(j)-interval/2,0,sigma);
    end
    transition(j,M)=1-normcdf(state(M)-rho*state(j)-interval/2,0,sigma);
end

%Check that each row of Transition Matrix sums to 1
if sum(transition')~=ones(1,M)
    nonzeros = find(sum(transition')-ones(1,M));
    disp('Error in Transition matrix');
    disp('the following rows do not sum to one:')
    disp(num2str(nonzeros));
end

%Calculate invariant distribution of Markov Chain
dist = ones(1,M)*(1/M);     %arbitrary initial distribution
diff = 1;
%Convergence of Distribution under sup norm
while diff>10^(-8)
    newdist=dist*transition;
    diff=max(abs(newdist-dist));
    dist=newdist;
end
invdist=dist;
    