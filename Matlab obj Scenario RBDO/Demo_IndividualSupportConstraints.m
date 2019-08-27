 % This script attempts computation of the support sets for indiviual
 % requirments and plot them for Case 1 Case 2 and Case 3
clc; clear 
close all 
%% Select Case Study
N=1e3; % number of available scenarios  
for Case_LSF=[1,2,3] % solve Case 1 2 and 3
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
%d*=SP3(alpha=0) multi performance, minimizes g_j \forall j=1,..,ng
Ng=size(Reliability_Dnominal.G,2);
Theta0sp3=[dn]; % set iniial condition   
alphasp3=0;
[Theta_opt3, Gamma_min_sp3_0, exitflag,output] = RBDO.Optimize_SP3(alphasp3, Theta0sp3);
Rel_Dopt_sp3=RBDO.Compute_ReliabilityMetrics(Theta_opt3);

 
%% Assess Robustness SP3, of individual requirements
Tol=1e-6; % tollerance of the new desing w.r.t. the reference optimal 
Deltaindex=[2,1];
if Case_LSF==1; PlotSub=1; elseif Case_LSF==2; PlotSub=3; elseif Case_LSF==3; PlotSub=5;end
subtightplot(3,2,PlotSub,0.01,0.1,0.05)
RBDO.plot_2D_SafeFailDomains_and_Scenarios(Theta_opt3,Deltaindex);
subtightplot(3,2, PlotSub+1,0.01,0.1,0.05);
RBDO.plot_2D_SafeFailDomains_and_Scenarios(Theta_opt3 ,Deltaindex); 

STYLE={'dg','*g' ,'og' ,'+g' };
STYLE2={'dr','*k' ,'*r' ,'dk' };

for G2examin=1:Ng
%Theta0=Theta_opt;
[ScenarioRobustness_add_SP3]=RBDO.ScenarioConstraintsSp3_addMethod(G2examin,Tol,alphasp3,Theta0sp3);
display(ScenarioRobustness_add_SP3)
Beta=1e-8;
epsilon=RBDO.getEpsilon(ScenarioRobustness_add_SP3.Cardinality,Beta);
 
subtightplot(3,2,PlotSub,0.01,0.1,0.05);
hold on
plot3(delta(ScenarioRobustness_add_SP3.Support_Set_Min,Deltaindex(1)),delta(ScenarioRobustness_add_SP3.Support_Set_Min,Deltaindex(2)),ones(ScenarioRobustness_add_SP3.Cardinality,1),STYLE{G2examin})
subtightplot(3,2, PlotSub+1,0.01,0.1,0.05);
plot3(delta(Rel_Dopt_sp3.G(:,G2examin)>0,Deltaindex(1)),...
      delta( Rel_Dopt_sp3.G(:,G2examin)>0,Deltaindex(2)),...
      ones(sum(Rel_Dopt_sp3.G(:,G2examin)>0),1),STYLE2{G2examin})

% save
IndividualRobustness{Case_LSF,G2examin}=ScenarioRobustness_add_SP3;
end
 
end 