%% Script-3: Run Scenario Program for different Rho Values
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',50000); % spq is needed to guaratnee high accuracy

%1) compute reliability for baseline desing
Rel_nominal=ComputeReliabilityPerformance(dn, delta,g_fun);  % get reliability

%2) run the scenario program optimizer for different rho values
Nrho=10;
VaR=0; % this defines the value at risk to compute the CVAR
RHOLINSPACE=linspace(0.0001,10,Nrho);

for i=1:Nrho
    
    RHO= RHOLINSPACE(i); % weight given to the violation
    [Results,X_opt,J_opt,exitflag,~]=ScenarioOptimizerCVaR(VaR,w_fun,delta,J_fun,RHO,dn,LBd,UBd,options);
    
    % get reliability performance for the delta_testing
    Rel_sp_opt=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);
    FailureProbability(i)=Rel_sp_opt.Pf_all;
    Pf_ind(:,i)=Rel_sp_opt.Pf_individual;
    Wmax(i)=Rel_sp_opt.w_max;
    Cost(i)=J_opt;
    Design(:,i)=Results.dopt;
    sN(i)=Results.Support.Size; % get size of the support sceanarios
    
    % Robustness
    beta=10^-8;
    OutWeJ = getWaitandJudgeEpsilon(sN(i),N ,beta);% get Wait-and-judge reliability
    epsilon_wait_and_judge =OutWeJ(end); % get Wait-and-judge reliability for the computed support size
    [EpsilonLU(1,i), EpsilonLU(2,i)] = epsLU(sN(i),N, beta);
    Epsi =Get_Epsilon_Individual_Requirements(g_fun,delta,Results.dopt,beta,VaR);
    Tempeps=Epsi';
    Epsilon_ind(:,i)=Tempeps(:);
end


 Plot_Bounds_and_groundtruth(EpsilonLU,RHOLINSPACE,FailureProbability)

for i=1:4
subplot(2,2,i);Plot_Bounds_and_groundtruth(Epsilon_ind([2*(i-1)+1,2*(i-1)+2],:),RHOLINSPACE,Pf_ind(i,:))
end
 