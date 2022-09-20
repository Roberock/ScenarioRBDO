clc; clear variables; close all
%% Example of truss strure RBDO
addpath([pwd '\G_functions'])
addpath([pwd '\G_functions\Truss'])
N=100; Ntest=10^4;
%% Example
% dn=[0.3075 -0.3164 -0.0973 -0.0188 -0.3164 0.5822 -0.0703 -0.0993 ...
%     -0.0973 -0.0703 0.2277 0.2661 -0.0188 -0.0993 0.2661 0.7100 -0.0191 ...
%     0.2733 -0.0920 0.4325 0.0803 -0.3821 0.4496 -0.2032]; % nominal design from previous studies
dn=[0.3075 -0.3164 -0.0973 -0.0188 0.5822 -0.0703 -0.0993  0.2277 0.2661  0.7100... % P matrix components
    -0.0191 0.2733 -0.0920 0.4325 0.0803 -0.3821 0.4496 -0.2032];% W matrix components
Nd=length(dn);
Ng=1;
LBd=dn-abs(dn)*0.3; UBd=dn+abs(dn)*0.3;

delta_mu=[ -2.93 -4.75 0.78 0.086 -0.11 0.1 -0.042 2.601 -0.29 -3.91 0.035 -2.5335 0.31];
DGM= @(N) normrnd(0,0.2 ,[N,13]).*delta_mu +delta_mu;


g_fun = @(d,delta)    G_fun_lateralmotioncontroller(d,delta);  % max eigneval
J_fun = @(d)          trace(reshape(d(1:16),4,4));  % Trace of P matrix
w_fun = @(d,delta)    max(g_fun(d,delta),[],2); % Worst-Case-Performance function

% data
delta=DGM(N);
delta_testing=DGM(Ntest);
%1) compute nomilan desing reliability
Rel_nominal_controller=ComputeReliabilityPerformance(dn, delta_testing,g_fun);  % get reliability

%2) run scenario program optimizer
VaR=0;
RHO=10^5; tic
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',25000,'display','iter','UseParallel',false);
[Results_controller,X_opt,~,~,optresults]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dn,LBd,UBd,options);
CompTime=toc;
Rel_sp_opt=ComputeReliabilityPerformance(Results_controller.dopt ,delta_testing,g_fun);  % get reliability
sN=Results_controller.Support.Size; % get size of the support sceanarios

% Robustness
beta=10^-4;
OutWeJ = getWaitandJudgeEpsilon(sN ,N,beta);% get Wait-and-judge reliability
epsilon_wait_and_judge =OutWeJ(end); % get Wait-and-judge reliability for the computed support size
[EpsilonLU(1), EpsilonLU(2)]= epsLU(sN,N, beta);
Epsilon_ind=Get_Epsilon_Individual_Requirements(g_fun,delta,Results_controller.dopt,beta,VaR);

disp(['True Failure Probability for the optimized design: Pf = '...
        num2str(Rel_sp_opt.Pf_all)])
   disp(['Number of support constraints in the scenario RBDO problem '...
        num2str(Results_controller.Support.Size)]) 
 disp(['Scenario bounds on the failure probability ['...
        num2str(EpsilonLU) ']'])  
 disp(['Obtained For a confidence level beta = '...
        num2str(beta) ' and number of scenario constraints N = ' num2str(N) ])      
%% Start with GA
% Nd=length(dn); % get number of desing variables
% Ndelta=size(delta,1);
% LB=[LBd, VaR*ones(1,Ndelta)]; % lower bounds on d and the slack variables \zeta_i
% UB=[UBd, inf*ones(1,Ndelta)]; %upper bounds on d and the slack variables \zeta_i
% X0=[dn, LB(Nd+1:end)]; % an initial guess
% obj= @(x) J_fun(x(1:Nd))+ RHO*sum(x(Nd+1:end)+LB(Nd+1:end)); % The cost function with slack variables J(d)+sum_i \zeta_i
% nonlincon= @(x) Scenario_Pf_constraint(x,delta,w_fun,Nd) ;
% %% run ga optimzer
% % Scenario_Pf_constraint compute the constraints defined as w(d,delta_i)>\zeta_i
% %% run fmincon optimzer
% options = optimoptions('ga','MaxTime',3600*1,'display','iter','UseParallel',true);
%
% Scenario_Pf_constraint compute the constraints defined as w(d,delta_i)>\zeta_i
% [X_opt,J_opt,exitflag,output]= ga(@(x) obj(x),... % objective function
%     length(X0),[],[],[],[], LB, UB,...% initial guess and design bounds
%     @(x) nonlincon(x),options); % non linear constraints and options

% Rel_sp_opt=ComputeReliabilityPerformance(X_opt(1:length(dn)) ,delta_testing,g_fun);  % get reliability


%% Example of 25 elements truss strur
% D=Data25; % Define the truss (geometry, connectivity, etc)
% LBd=D.LB; UBd=D.UB;
% dn= [0.3, 2.00, 3.40, 0.1, 0.1 , 1.20, 1.9, 3.40]; % Define a random design (section area for each group of elements)
% Nd=length(dn);
% Mu_load= D.LoadCase_rnd(D.LoadCase_rnd~=0);
% DGM= @(N)  [Mu_load'+abs(Mu_load)'.*trnd(4,[N,7])/50, normrnd(0,10^6,[N,25])];
% g_fun = @(d,delta) g_fun_ST25C(d,delta,D);
% J_fun = @(d,delta) J_fun_ST25C(d,D);
% w_fun = @(d,delta) max(g_fun(d,delta),[],2); % Worst-Case-Performance function
%
% % data
% delta=DGM(N);
% delta_testing=DGM(Ntest);
% %1) compute nomilan desing reliability
% Rel_nominal=ComputeReliabilityPerformance(dn, delta,g_fun);  % get reliability
%
% %2) run scenario program optimizer
% VaR=0; RHO=10^4; tic
% options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',15000,'display','iter','UseParallel',true);
% [Results,X_opt,J_opt,~,optresults]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dn,LBd,UBd,options);
% CompTime=toc;
% Rel_sp_opt=ComputeReliabilityPerformance(Results.dopt ,delta,g_fun);  % get reliability
% sN=Results.Support.Size; % get size of the support sceanarios
%
% % Robustness
% beta=10^-8;
% OutWeJ = getWaitandJudgeEpsilon(sN ,N,beta);% get Wait-and-judge reliability
% epsilon_wait_and_judge =OutWeJ(end); % get Wait-and-judge reliability for the computed support size
% [EpsilonLU(1), EpsilonLU(2)]= epsLU(sN,N, beta);
% Epsilon_ind=Get_Epsilon_Individual_Requirements(g_fun,delta,Results.dopt,beta,VaR);

