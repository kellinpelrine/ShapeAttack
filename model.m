function returntest = model(gamma,epsilon,alpha,rho,sigmaerror,mu,IYssratio,KYssratio,nss,kss,delta,beta,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rbc_pset.m
%
% Calibration and Evaluation of the RBC Model
% Problem Set 1
% written by YOUR NAME
%
% calls functions u.m and markovchain.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tic
%clear
%%global beta gamma b epsilon delta alpha M z ztrans
% 
% %% Set Parameter Values
% gamma       = 1;                %assuming log utility
% epsilon     = 2;                %elasticity of labor supply
% alpha       = .4;               %capital's share of income
% rho         = .95;
% sigmaerror  =.01;
% mu          = 0;                %long-run growth rate
% %Steady State Values
% IYssratio   = .25;
% KYssratio   = 14;
% nss         = 1/3;
% kss         = 14^(.06) / 3;              
% %Calibrated Parameter Values
% delta   = 1/56;                  %depreciation rate
% beta    = .964;                  %discount rate
% b       = 7.2;                  %disutility of labor multiple


%% Create Grids
% Grid for Productivity
% productivity evolves according to an M-state Markov Process
% log(z_t) = rho*log(z_(t-1)) + error_t+1
M         = 12;
suppwidth = 3;
[logzstate, ztrans, invzdist]=markovchain(M, rho, sigmaerror, suppwidth);
z=exp(logzstate);

% Grid for Capital
khigh   =kss*1.2;
klow    =kss*0.8;
kap     =[klow:.05:khigh]';
Nkap    =length(kap);

% Grid for labor
labor   =[0.2:.01:0.7];
Nlab    =length(labor);

% Declare Matrices
v       = zeros(Nkap,M);    %initial guess for v
kpolicy = zeros(Nkap,M);
Utilde  = zeros(Nkap,Nkap,M);
npolicy = zeros(Nkap,Nkap,M);

%% Create Consumption and Utility Matrices
% first dimension: tomorrow's k     (index i)
% second dimension: today's k       (index j)
% third dimension: productivity     (index m)
% fourth dimension: labor
for i=1:Nkap
    for j=1:Nkap
        for m=1:M
            cons=(1-delta)*kap(j)+z(m)*(kap(j)^alpha)*(labor.^(1-alpha))-kap(i);
            C=max(cons,.0001);
            Utility=u(C,labor,gamma,b,epsilon);

            s = find(cons<=0);              % infeasible consumption choice
            Utility(s) = -1e20;

            % maximize over labor, reduce U to 3-dimensional array Utilde
            [maxu,nchoice]= max(Utility);
            % maxv contains the maximum value in this vector
            % nchoice contains the index of the maximum value
            Utilde(i,j,m)     =maxu;
            npolicy(i,j,m)    =nchoice;
        end
    end
end

%disp('start iteration')
%% Value Function Iteration
%initialize iteration
crit = 1; iter = 0;
vnew = zeros(Nkap,M);
kpolicynew = zeros(Nkap,M);
npolicynew = zeros(Nkap,M);
%iterations:
while crit~=0
    iter=iter+1;
    
    for m=1:M
        % Expected Value tomorrow for each capital chosen
        futureval = ztrans(m,:)*(v');  %1xNkap vector: each entry is EV(k)
        
        [maxv,kchoice]= max(Utilde(:,:,m)+beta*futureval'*ones(1,Nkap));
        % maxv is a 1xNkap row vector containing the maximum element from each
        % column of Utilde(:,:,m)+beta*futureval'*ones(1,Nkap)
        % kchoice is a 1xNkap row vector containing the indices of the
        % maximum values of Utilde(:,:,m)+beta*futureval'*ones(1,Nkap)
        
        vnew(:,m)=maxv';
        kpolicynew(:,m)=kchoice';
    end
    
    % criterion: policy function convergence
    critnew=max(any(kpolicynew-kpolicy));  %returns 0 if policy function is same as before
    v=vnew;
    kpolicy=kpolicynew;
    crit=critnew;
end

for m=1:M
    for j=1:Nkap
        npolicynew(j,m)=npolicy(kpolicy(j,m),j,m);
    end
end
% Write policy function in terms of capital
gcapital    =kap(kpolicy);
glabor      =labor(npolicynew);

gcapsmooth      =spline(z,gcapital);
glabsmooth      =spline(z,glabor);

%% Monte-Carlo
shock       =sigmaerror*randn(100,1000);
prodSim     =zeros(100,1000);
kSim        =zeros(100,1000);
nSim        =zeros(100,1000);
qSim        =zeros(100,1000);
%cSim        =zeros(100,500);
iSim        =zeros(100,1000);
%wSim        =zeros(100,500);
%rSim        =zeros(100,500);
%stdSim      =zeros(2,500);
corrSim     =zeros(2,2,1000);
for t=1:100
    if t==1
        prodSim(1,:)=exp(shock(1,:));
    else
        prodSim(t,:) =exp(rho*log(prodSim(t-1,:))+shock(t,:));
    end
end

for s=1:1000
    kSimTemp = zeros(100,1);
    iSimTemp = zeros(100,1);
    nSimTemp = zeros(100,1);
    qSimTemp = zeros(100,1);
    cSimTemp = zeros(100,1);
    for t=1:100
        if t==1
            kSimTemp(t) = kss;
            iSimTemp(t) = delta*kss;
            %kSim(t,s)  =kss;
            %iSim(t,s)  =delta*kss;
        else
            kvec       =ppval(gcapsmooth,prodSim(t,s));
            kSimTemp(t) =spline(kap,kvec,kSimTemp(t-1));
            iSimTemp(t) =kSimTemp(t)-(1-delta)*kSimTemp(t-1);
            %kSim(t,s)  =spline(kap,kvec,kSim(t-1,s));
            %iSim(t,s)  =kSim(t,s)-(1-delta)*kSim(t-1,s);
%            if kSimTemp(t) < 0
%                kSimTemp(t) = 0
%            end
%             if any(isinf(kSimTemp(t)),'all')
%                 kSimTemp(1:t-1)
%                 kvec
%                 return
%             end
        end
        nvec       =ppval(glabsmooth,prodSim(t,s));
        nSimTemp(t)  =spline(kap,nvec,kSimTemp(t));
        qSimTemp(t)  =prodSim(t,s)*(kSimTemp(t)^alpha)*(nSimTemp(t)^(1-alpha));
        cSimTemp(t)  =qSimTemp(t)-iSimTemp(t); 
        
        %nSim(t,s)  =spline(kap,nvec,kSim(t,s));
        %qSim(t,s)  =prodSim(t,s)*(kSim(t,s)^alpha)*(nSim(t,s)^(1-alpha));
        %cSim(t,s)  =qSim(t,s)-iSim(t,s);
        %wSim(t,s)  =(1-alpha)*(qSim(t,s)/nSim(t,s));
        %rSim(t,s)  =alpha*(qSim(t,s)/kSim(t,s));
    end
%     kSim(:,s) = kSimTemp;
%     iSim(:,s) = iSimTemp;
%     nSim(:,s) = nSimTemp;
%     qSim(:,s) = qSimTemp;
%     cSim(:,s) = cSimTemp;
    % calculate correlations and standard deviations
    %M              =[qSim(:,s) nSim(:,s) kSim(:,s) cSim(:,s) iSim(:,s) wSim(:,s) rSim(:,s)];
    M = [qSimTemp cSimTemp];    
    corrSim(:,:,s) =corr(M);
    %stdSim(:,s)    =std(M)';
end

%Standard Deviations and Correlations
Corr = nanmean(corrSim,3);
%Corr    =sum(corrSim,3)/500;
%Stdev   =[sum(stdSim,2)/500]';
returntest = real(Corr(1,2));
 
% %Plot an artifical sequence
% time            = [1:100]'; seq=14; %arbitrary sequence
% figure(1);
% subplot(5,1,1); plot(time,prodSim(:,seq)); xlabel('t'); ylabel('z');
% subplot(5,1,2); xxx %add plot for k here
% subplot(5,1,3); xxx %add plot for c here
% subplot(5,1,4); xxx %add plot for q(=y) here
% subplot(5,1,5); xxx %add plot for n here
% 
% clear time shock prodSim kSim nSim qSim cSim iSim wSim rSim stdSim corrSim
% 
% 
% %% Generate and Plot Impulse Responses
% nperiods    =60;
% prod        =zeros(1,nperiods);
% kIR         =zeros(1,nperiods);
% nIR         =zeros(1,nperiods);
% qIR         =zeros(1,nperiods);
% iIR         =zeros(1,nperiods);
% cIR         =zeros(1,nperiods);
% wIR         =zeros(1,nperiods);
% rIR         =zeros(1,nperiods);
% 
% shock           =sigmaerror; %single deviation shock
% for t=1:nperiods
%     if t==1
%         prod(t) =exp(shock);
%         kIR(t)  =kss;
%         iIR(t)  =delta*kss;
%     else
%         prod(t) =exp(rho*log(prod(t-1)));
%         kvec    =ppval(gcapsmooth,prod(t-1));
%         kIR(t)  =spline(kap,kvec,kIR(t-1));
%         iIR(t)  =kIR(t)-(1-delta)*kIR(t-1);
%     end
%     nvec    =ppval(glabsmooth,prod(t));
%     nIR(t)  =spline(kap,nvec,kIR(t));
%     qIR(t)  =prod(t)*(kIR(t)^alpha)*(nIR(t)^(1-alpha));
%     cIR(t)  =qIR(t)-iIR(t);
%     wIR(t)  =(1-alpha)*(qIR(t)/nIR(t));
%     rIR(t)  =alpha*(qIR(t)/kIR(t));
%     
% end
% 
% time            = [0:nperiods-1]';
% figure(2); xxx %plot Productivity Shock 
% figure(3); xxx %plot Impulse Response of Output
% figure(4); xxx %plot Impulse Response of Employment
% figure(5); xxx %plot Impulse Response of Consumption
% figure(6); xxx %plot Impulse Response of Investment
% figure(7); xxx %plot Impulse Response of Capital
% figure(8); xxx %plot Impulse Response of Wage Rate
% figure(9); xxx %plot Impulse Response of Interest Rate

%toc