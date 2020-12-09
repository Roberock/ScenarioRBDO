%% Script-1: Run Scenario Program optimizer and compare solution ot the nominal case and to alternative CVaR minimization algorithm
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',50000,'Display','iter'); % spq is needed to guaratnee high accuracy
RHO= 100; % weight given to the violation
VaR=0; % this defines the value at risk to compute the CVAR

%1) compute nomilan desing reliability
Rel_nominal=ComputeReliabilityPerformance(dn, delta,g_fun);  % get reliability

%1.1) Chance-Constrained CVAR problem
alpha=0.99;
[Resultscvar,X_opt_cvar,J_opt_CVAR,~,~]=CVaR_CCP(delta,alpha,w_fun,J_fun,dn,LBd,UBd);
Rel_CVaR_opt=ComputeReliabilityPerformance(Resultscvar.dopt ,delta,g_fun);  % get reliability

%2) run scenario program optimizer
[Results,X_opt,J_opt,~,~]=ScenarioOptimizerCVaR(VaR,w_fun,delta,J_fun,RHO,dn,LBd,UBd,options);
Rel_sp_opt=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);  % get reliability
sN=Results.Support.Size; % get size of the support sceanarios

% Robustness
beta=10^-8;
OutWeJ = getWaitandJudgeEpsilon(sN ,N,beta);% get Wait-and-judge reliability
epsilon_wait_and_judge =OutWeJ(end); % get Wait-and-judge reliability for the computed support size
[EpsilonLU(1), EpsilonLU(2)]= epsLU(sN,N, beta);
Epsilon_ind=Get_Epsilon_Individual_Requirements(g_fun,delta,Results.dopt,beta,VaR);
% DisplayTableResults
 