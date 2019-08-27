% TUTORIAL: ScenarioRBDO class
% The ScenarioRBDO Object is a collection of methods to perform Scenario Reliability-Based-Design-Optimization
  
clc
clear 
close all

%% Select Case Study
N=1e3; % number of available scenarios 
Case_LSF=8; % select the case study
%Case_LSF=1  g in R^2 linear competitive LSF, delta in R^2
%Case_LSF=2  g in R^2 non-linear convex/concave performance functions Groteman et al. 2011, delta in R^2
%Case_LSF=3  g in R^4 non-linear performance functions paralel system, delta in R^2
%Case_LSF=4  g in R^6 Speed Reducer modified, delta in R^6, d in R^5
%Case_LSF=5  g in R^6 Antenna tower reliability problem, delta in R^27
%Case_LSF=6  g in R^3 Robust Control Benchmark, delta in R^2 % available N=1e3;
%Case_LSF=7; g in R^33, delta in R^34, N=1e4 Ndes=user-defined, is a modified Rosenbrock function (optimization benchmark with many local min/max) 
%Case_LSF=8; g in R^2 sinusoidal non linear 

[g_fun,delta,dn,LBd,UBd]=Select_Case_Study(N,Case_LSF);

%%  Create an object SCENARIO RBDO OBJECT
OptimizerData.LB=[LBd];
OptimizerData.UB=[UBd];
OptimizerData.options= optimoptions('fmincon','Display','off','OptimalityTolerance',1e-6);
[OptimizerData.A, OptimizerData.B, OptimizerData.Aeq, OptimizerData.Beq]=deal([]);
RBDO=ScenarioRBDO('delta',delta,'dn',dn,'g_fun',g_fun,'OptimizerData',OptimizerData);

%% Anlyse Reliability of the nominal desing
Reliability_Dnominal=RBDO.Compute_ReliabilityMetrics(dn);
display(Reliability_Dnominal)
%RBDO.plot_2D_SafeFailDomains_and_Scenarios(dn,[1,2],50) % plot failure regions for g1 g2

%% RUN Scenario Optimizations and analyse d* performance:
% d*=SP1(alpha=0) minimize worst case w=max_j g_j 
alpha=0; % define 1-\alpha percentile to keep \alpha=0 means all the scenarios
Theta0=[dn];
% Gamma_min not equal w max as individual requirements are normalized in [-1,+1]
[Theta_opt, Gamma_min_sp1_0, exitflag_sp1_0,output_sp1_0] = RBDO.Optimize_SP1(0, Theta0);
Rel_Dopt=RBDO.Compute_ReliabilityMetrics(Theta_opt);

% d*=SP1(alpha=0.05) minimize worst case percentile: p_95[w]=F^{-1}_w(1-\alpha)
[Theta_opt05, Gamma_min_sp1_05, ~,~] = RBDO.Optimize_SP1(0.05, Theta0);
Rel_Dopt05=RBDO.Compute_ReliabilityMetrics(Theta_opt05);

% d*=SP2() minimize given data failure probability
[Theta_optpf,Pf_min] = RBDO.Optimize_SP2();
Rel_Doptpf=RBDO.Compute_ReliabilityMetrics(Theta_optpf);

%d*=SP3(alpha=0) multi performance, minimizes g_j \forall j=1,..,ng
Ng=size(Reliability_Dnominal.G,2);
Theta0sp3=[dn]; % set iniial condition   
alphasp3=0;
[Theta_opt3, Gamma_min_sp3_0, exitflag,output] = RBDO.Optimize_SP3(alphasp3, Theta0sp3);
Rel_Dopt_sp3=RBDO.Compute_ReliabilityMetrics(Theta_opt3);


%% Plots, some examples:
% compare failure safe regions, nominal vs optimal desings
Deltaindex=[2,1];
subtightplot(1,2,1,0.01,0.1,0.05)
RBDO.plot_2D_SafeFailDomains_and_Scenarios(dn,Deltaindex);
subtightplot(1,2,2,0.01,0.1,0.05);
RBDO.plot_2D_SafeFailDomains_and_Scenarios(Theta_opt,Deltaindex);
  % show individual failure regions for each requirement
  Gindex=[1,2];
  RBDO.plot_2D_SafeFailDomains_and_Scenarios(dn,Deltaindex,Gindex)
  RBDO.plot_2D_SafeFailDomains_and_Scenarios(Theta_opt,Deltaindex,Gindex)
%% compare show g(z-axes) vs delta, nominal vs optimal desings
Gindex=[1,2];
RBDO.plot_deta_vs_G(Theta_opt,Deltaindex,Gindex);
% RBDO.plot_deta_vs_G(Theta_opt,Deltaindex,Gindex)
%% show gi vs delta indexed, INPUTS (d,idx)  d=desing idx= index of g_i to be sorted (for better graphical output)
figure
G2Sort=1; % index of the requirement to be sorted
RBDO.plot_detaIndex_vs_G(Theta_opt,G2Sort) %
 

%% ROBUSTNESS ASSESSMENT METHODS % (Can Be Time Consuming)

%% Assess Robustness SP1
Tol=1e-6;
%Theta0=Theta_opt;
[ScenarioRobustness_add]=RBDO.ScenarioConstraints_addMethod(Tol,alpha,Theta0);
display(ScenarioRobustness_add)
Beta=1e-8;
epsilon=RBDO.getEpsilon(ScenarioRobustness_add.Cardinality,Beta);

%% Assess Robustness SP3, of individual requirements
Tol=1e-6; % tollerance of the new desing w.r.t. the reference optimal 
G2examin=2; % it examins requirement g_2
[ScenarioRobustness_add_SP3]=RBDO.ScenarioConstraintsSp3_addMethod(G2examin,Tol,alpha,Theta0sp3);
display(ScenarioRobustness_add_SP3)
 
%% Assess Robustness method 2 
%[ScenarioRobustness_remove]=RBDO.ScenarioConstraints_removeMethod(Tol,alpha,Theta0); 
%% Assess Robustness for program SP3, individual requirements
% G2examin=1;
% [Rbst]=RBDO.ScenarioConstraintsSp3_removeMethod(G2examin,Tol,alphasp3,Theta0sp3);
%% Plots, some example: 
%% plot in the d-w plane the scenario constraints INPUTS (d,idx) d=desing idx= index of delta to plot
% Scenario2Show=[1:1:100]; % show first 5 scenario constraints
Scenario2Show=ScenarioRobustness_add.Support_Set_Min; % to check only the supports 
RBDO.Plot_ScenarioConstraints(Theta_opt,Scenario2Show)
 
%% Improve Design: Remove Scenarios
Support_Set=ScenarioRobustness_add.Support_Set_Min;
%RemovalType=0;% remove one at a time from the set of support S with replacement
%RemovalType =1 % remove one at a time from the set of support S without replacement (an order of removal should be provided)
RemovalType=2;% remove all together
[ScenarioRobustness_removeS]=RBDO.RemoveConstraints(Support_Set,ScenarioRobustness_add,RemovalType);
Rel_Dopt_removeS=RBDO.Compute_ReliabilityMetrics(ScenarioRobustness_removeS.Theta_opt_i);
