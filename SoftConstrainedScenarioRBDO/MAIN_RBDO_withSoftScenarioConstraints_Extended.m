clc; clear variables; close all
%% CONSIDER THE FOLLOWING SCENARIO PRORAMs WITH SOFT-CONSTRAINTS: 
% Program-1
% min_{d\in \Theta , \zeta^{(i)}>0} \ lbrace J(d) +\rho \sum\limits_{i=1}^{N}  \zeta^{(i)}
% Such that: w(d,\delta{(i)} \leq \zeta^{(i)} \rbrace
% where $\rho$ is a parameter weighting the cost of violation for w and
% w(d,\delta)=\max\limits_{j\in\{1,..,n_g \}| g_j(d,\delta) %worst-case reliability performance function
% Program-2
% min_{d\in \Theta , \zeta_j^{(i)}>0} \ lbrace J(d) +\sum\limits_{j=1}^{n_g} \rho_j \sum\limits_{i=1}^{N}  \zeta_j^{(i)}
% Such that: g_j(d,\delta) \leq \zeta_j^{(i)} i=1,...,N,~j=1,..,n_g\rbrace
% where $\rho_j$ are parameters weighting the cost of violation on the reliability requirement g_j

% for this probelm the support scenarios (complexity0 will be equal to
% S_N^*= the number of active constriaints + the number of violating constraints
%% START with PROGRAM-1
%% Load Case study
CASESTUDY=2;
N=200;  % available samples
[g_fun,delta,dn,LBd,UBd,DGM,Nd,Ng]=Select_Convex_Case_Study(N,CASESTUDY);
N_testing=10^5;  % number of scenarios for validation
delta_testing=DGM(N_testing); % generate more scenarios for validation
J_fun= @(d) sum(d); % define cost function
% g_fun=@(d,delta) g_fun(d,delta)./(1+abs(g_fun(d,delta)))
w_fun=@(d,delta) max(g_fun(d,delta),[],2); % Worst-Case-Performance function
%w_fun=@(d,delta) w_fun(d,delta)./(1+abs(w_fun(d,delta)));
%% Script-1: Run Scenario Program optimizer and compare solution ot the nominal case and to alternative CVaR minimization algorithm
%% Experiment 1: SP run for for \lambda=0 and \rho=100 compare to nominal desing
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',50000);
Script1
%Table_1
%% Script-2:  run for increasing Number of samples, N=linspace(,,)
Script2
%Table_2
%% Script-3:  run for increasing Number of samples, N=linspace(,,)
Script3
%% Script 4: OPTIMIZE for N=linspace() lambda=linspace( ) and rho=100
Script4
 
%% Experiment 3: OPTIMIZE for lambda=0 and rho=linspace(0,100,50)
Nrhos=50;
RHOvals= linspace(0,100,Nrhos) ;
[AlphaValueAtRisk,CVAR,Pf_all,Jopt,SupportSize,epsilon_rho_convex,epsilon_rho_nnconvex,Alpha_true,CVAR_true,Pf_true]=deal(zeros(1,Nrhos));
for i=1:length(RHOvals)
    display(num2str(i));
    VaR=0; % this defines the value at risk to compute the CVAR
    [Results,X_opt,J_opt,exitflag,output]=ScenarioOptimizerCVaR(VaR,w_fun,delta,J_fun,RHOvals(i),dn,LBd,UBd,options);
    % examin desing reliability
    Rel=ComputeReliabilityPerformance(Results.dopt ,delta,g_fun);  % get reliability
    AlphaValueAtRisk(i)=Results.AlphaValueAtRisk;
    dopt_exp3(i,:)=Results.dopt;
    CVAR(i)=CVaR(Rel.w,1-Rel.Pf_all);
    Pf_all(i)=Rel.Pf_all;
    Jopt(i)=Results.Jopt;
    SupportSize(i)=Results.Support.Size;
    
    % examin scenario prescpective reliability
    OutWeJ =getWaitandJudgeEpsilon(SupportSize(i),N,1e-8);
    epsilon_rho_convex(i)=OutWeJ(end);
    epsilon_rho_nnconvex(i) =getConfidence_nonconvex(SupportSize(i),1e-8,N);
    
    %  true failure probability
    res1e5=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);
    Alpha_true(i)=mean(res1e5.w>VaR);
    CVAR_true(i)=mean(res1e5.w(res1e5.w>=VaR));
    Pf_true(i)=mean(res1e5.w>0);
end

figure()
yyaxis left
%plot(AlphaValueAtRisk,epsilon_rho_nnconvex)
hold on
plot(RHOvals,epsilon_rho_nnconvex,'LineWidth',2)
plot(RHOvals,Pf_all,'LineWidth',2)
grid on;box on
xlabel('\rho ');ylabel(' \epsilon(s_N^*); P_f(d^*)')

yyaxis right
plot(RHOvals, Jopt)
%plot(AlphaValueAtRisk, LambdaVal,'LineWidth',1)
plot(RHOvals,CVAR,'LineWidth',1)
grid on;box on
ylabel('VaR ; CVaR ') ; legend('\epsilon(s_N^*)','P_f(d^*)','VaR','CVaR(d^*)')
set(gca,'FontSize',14)

%% Analyse individual requirements

%% Program-2
% min_{d\in \Theta , \zeta_j^{(i)}>0} \ lbrace J(d) +\sum\limits_{j=1}^{n_g} \rho_j \sum\limits_{i=1}^{N}  \zeta_j^{(i)}
% Such that: g_j(d,\delta) \leq \zeta_j^{(i)} i=1,...,N,~j=1,..,n_g\rbrace
% where $\rho_j$ are parameters weighting the cost of violation on the reliability requirement g_j

%CASESTUDY=4;
N=500;  % available samples
[g_fun,delta,dn,LBd,UBd,DGM,Nd,Ng]=Select_Convex_Case_Study(N,CASESTUDY);
Ndelta=size(delta,1);
LB=[LBd]; UB=[UBd];
VaR=0;
for i=1:Ng
    LB=[LB, VaR*ones(1,Ndelta)]; % lower bounds on d and the slack variables \zeta_i
    UB=[UB, inf*ones(1,Ndelta)]; %upper bounds on d and the slack variables \zeta_i
end
X0=[dn, LB(Nd+1:end)]; % an initial guess

for rj=1:10
    RHO_j=[1, 1, rj/5, 1];
    obj= @(x) J_fun(x(1:Nd));
    for i=1:Ng
        obj=  @(x) obj(x)+  RHO_j(i)*sum(x(Nd+1+(i-1)*N:Nd+N+(i-1)*N)+LB(Nd+1+(i-1)*N:Nd+N+(i-1)*N));
    end
    
    %  run fmincon optimzer
    options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',50000,'Display','iter-detailed');
    
    % Scenario_Pf_constraint compute the constraints defined as w(d,delta_i)>\zeta_i
    g_norm =  @(d,delta) g_fun(d,delta)./max(g_fun(dn,delta));
    [X_opt,J_opt,exitflag,output]= fmincon(@(x) obj(x),... % objective function
        X0,[],[],[],[], LB, UB,...% initial guess and design bounds
        @(x) Scenario_Pf_constraint_multipleSoftG(x,delta,g_fun,Nd),options); % non linear constraints and options
    
    %  evaluate reliability
    Rel_multiRHOG=ComputeReliabilityPerformance(X_opt(1:Nd) ,delta_testing,g_fun);
    
    %  Collect Results
    Jopt(rj)=J_fun(X_opt(1:Nd));
    dopt(:,rj)= X_opt(1:Nd);
    Results.Zopt= reshape(X_opt(Nd+1:end),size(delta,1),Ng);
    Pf_all(rj)=Rel_multiRHOG.Pf_all;
    Pf_ind(:,rj)=Rel_multiRHOG.Pf_individual;
    Results.Support.Size=sum(g_fun(X_opt(1:Nd),delta)>=VaR);
    Results.Support.Scenarios= (g_fun(X_opt(1:Nd),delta)>=VaR);
    
    %robustness
    beta=10^-8;
    Epsi =Get_Epsilon_Individual_Requirements(g_fun,delta,dopt(:,rj),beta,0);
    Tempeps=Epsi';
    Epsilon_ind_multirho(:,rj)=Tempeps(:);
    
    
end

% plot
for i=1:Ng
    subplot(ceil(Ng/2),2,i);
    Plot_Bounds_and_groundtruth(Epsilon_ind_multirho([2*(i-1)+1,2*(i-1)+2],:),1:1:10,Pf_ind(i,:))
end



