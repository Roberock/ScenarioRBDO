%% DEMO: Solve problem 'Case 1' using the Class ScenarioRBDO
clc; clear variables;close all
 
%% define the reliability problem, available scenarios (delta) and performance function handle (g_fun)
dn=[1,1,1]; % assume nominal desing
N=5e2; % just 500 observations from the uncertainty space (scenarios)
delta=[normrnd(0,1,[N,1]),normrnd(0,2,[N,1])]; % normal distributed observations with mean 0 and standard dev [1,2] 
% g is defined such that g>0 is failure g<0 is safe g=0 is the limit state function;
g_fun1=@(delta,d) +1/d(1)*(delta(:,2))+1/d(2)*(delta(:,1))-d(3); % requirement 1
g_fun2=@(delta,d) +d(1)*(delta(:,1))-1/d(2)*(delta(:,2))-d(3); % requirement 2
g_fun=@(delta,d) [g_fun1(delta,d) g_fun2(delta,d)]; % the performance function 

%%  %  %  % Now define the optimization problem, bounds, constraints, and options %  %  %  %  %
OptimizerData.LB=[0.5,0.5,0.5]; %  lower bounds on the [d] search space
OptimizerData.UB=[2,2,2]; %  upper bounds  on the [d] search space
OptimizerData.options= optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-6);
[OptimizerData.A, OptimizerData.B, OptimizerData.Aeq, OptimizerData.Beq]=deal([]);

%%  %  %  %  %  Construct scenario reliability object %  %  %  %  %%  %  %  %  %
RBDO=ScenarioRBDO('delta',delta,'dn',dn,'g_fun',g_fun,'OptimizerData',OptimizerData);

%  %  % Assess reliability of the nominal desing %  %  %  %  %  %  %  %
Rel_dn=RBDO.Compute_ReliabilityMetrics(dn); display(Rel_dn) 

%%  %  %  %  %  RUN Optimization programs %  %  %  %  %%  %  %  %  %
%% SP1(alpha=0) %  % Optimize the desing by solving program SP1:%  %  %
% SP1(alpha=0) minimizes the worst case performance among all the scenarios %  %  %  %  %  %  %
alpha=0;  Theta0=[dn];
[Theta_opt_SP1_0, ~, ~, ~] = RBDO.Optimize_SP1(alpha, Theta0);
Rel_sp1_0=RBDO.Compute_ReliabilityMetrics(Theta_opt_SP1_0); display(Rel_dn) % evaluate the optimaldesing reliability

%% SP1(alpha=0.05)  %  % Optimize the desing by solving program SP1:%  %  %
% SP1(alpha=0.05) minimizes the worst case performance removing 5 % of the worst case scenarios % % % % % % %
alpha=0.05;  %gamma0=0; Theta0=[dn gamma0];
[Theta_opt_SP1_005, ~, ~, ~] = RBDO.Optimize_SP1(alpha, Theta0);
Rel_sp1_05=RBDO.Compute_ReliabilityMetrics(Theta_opt_SP1_005); display(Rel_sp1_05) % evaluate the optimaldesing reliability

%% SP2()  %  % Optimize the desing by solving program SP1:%  %  %
% SP2 minimizes the given-data estimator of the failure probability % % % % % % %
addpath([pwd '/GA_solver'])
[Theta_opt_SP2,Pf_min] = RBDO.Optimize_SP2();
Rel_sp2=RBDO.Compute_ReliabilityMetrics(Theta_opt_SP2); display(Rel_sp2) % evaluate the optimaldesing reliability

%% SP3(alpha=0)  %  % Optimize the desing by solving program SP3:%  %  %
% SP3(alpha=0) minimizes the norm of the worst cases for the ng requirement
addpath([pwd '/GA_solver'])
alpha=0;  %gamma0=zeros(1,RBDO.Ng); Theta0=[dn gamma0]; 
[Theta_opt_SP3_0,Pf_min] = RBDO.Optimize_SP3(alpha, Theta0);
Rel_sp3_0=RBDO.Compute_ReliabilityMetrics(Theta_opt_SP3_0); display(Rel_sp3_0) % evaluate the optimaldesing reliability
%% SP3(alpha=0.05)  %  % Optimize the desing by solving program SP3:%  %  %
alpha=0.05;  %Theta0SP3=[dn ];
[Theta_opt_SP3_05,Pf_min] = RBDO.Optimize_SP3(alpha, Theta0);
Rel_sp3_05=RBDO.Compute_ReliabilityMetrics(Theta_opt_SP3_05); display(Rel_sp3_05) % evaluate the optimaldesing reliability
%% SP3(alpha=0.2)  %  % Optimize the desing by solving program SP3:%  %  %
alpha=0.2;  %Theta0SP3=[dn ];
[Theta_opt_SP3_2,Pf_min] = RBDO.Optimize_SP3(alpha, Theta0);
Rel_sp3_2=RBDO.Compute_ReliabilityMetrics(Theta_opt_SP3_2); display(Rel_sp3_2) % evaluate the optimaldesing reliability


%% PLOT  SP3(0.05) vs SP3(0.2) vs SP3(0)
Deltaindex=[1,2];
subtightplot(1,3,1,0.01,0.1,0.05)
RBDO.plot_2D_SafeFailDomains_and_Scenarios(Theta_opt_SP3_0,Deltaindex);
subtightplot(1,3,2,0.01,0.1,0.05);
RBDO.plot_2D_SafeFailDomains_and_Scenarios(Theta_opt_SP3_05,Deltaindex);
subtightplot(1,3,3,0.01,0.1,0.05);
RBDO.plot_2D_SafeFailDomains_and_Scenarios(Theta_opt_SP3_2,Deltaindex);

%% Demo Robustness Analysis SP1(0)
%  %  % Robustnes of the optimal desing %  %  %  %  %
Tol=1e-6; % tollerance  the removed scenario
alpha=0;
[Robustness_optsp0]=RBDO.ScenarioConstraints_addMethod(Tol,alpha,Theta0); % find support scenarios using the add-method
display(Robustness_optsp0) 

%% Demo Robustness Analysis SP1(0.05)
%  %  % Robustnes of the optimal desing %  %  %  %  %
Tol=1e-6; % tollerance  the removed scenario
alpha=0.05;
[Robustness_optsp05]=RBDO.ScenarioConstraints_addMethod(Tol,alpha,Theta0); % find support scenarios using the add-method
display(Robustness_optsp05) 
  %%  %  % Robustnes of the optimal desing %  %  %  %  %
% Tol=1e-6; % tollerance  the removed scenario
% [Robustness_optsp0]=RBDO.ScenarioConstraints_removeMethod(Tol,alpha,Theta0); % find support scenarios using the add-method
