%% Script-2: Run Scenario Program for increasing number of samples N
RHO= 1; % weight given to the violation
VaR=0; % this defines the value at risk to compute the CVAR
NSAMPLES=[100:25:1000]; 
[CVAR95,CVAR_Fzero,Pf_all,Wmax,Jopt,SupportSize,epsilon_rho_convex]=deal(zeros(1,length(NSAMPLES)));
dopt_exp2=zeros(Nd,length(NSAMPLES));
[Gmax,Pf_ind]=deal(zeros(Ng,length(NSAMPLES)));
Epsilon_ind=zeros(Ng*2,NSAMPLES);
for i=1:length(NSAMPLES)
    delta_N=DGM(NSAMPLES(i));
    
    % run the optimizer
    [Results,X_opt,J_opt,exitflag,output]=ScenarioOptimizerCVaR(VaR,w_fun,delta_N,J_fun,RHO,dn,LBd,UBd,options);
     dopt_exp2(:,i)=Results.dopt;
     Jopt(1,i)=Results.Jopt;
     SupportSize(1,i)=Results.Support.Size;
      
    % get reliability  scores for the test set of scenarios
    Rel=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);
    CVAR_Fzero(1,i)=CVaR(Rel.w,1-Rel.Pf_all);
    CVAR95(1,i)=CVaR(Rel.w,0.95);
    Pf_all(1,i)=Rel.Pf_all;
    Pf_ind(:,i)=Rel.Pf_individual;
    Wmax(i)=Rel.w_max;
    Gmax(:,i)=Rel.g_max;

    % examin scenario prescpective reliability
    OutWeJ =getWaitandJudgeEpsilon(SupportSize(1,i),NSAMPLES(i),1e-8);
    epsilon_rho_convex(1,i)=OutWeJ(end);
    [EpsilonLU(1,i), EpsilonLU(2,i)]= epsLU(SupportSize(1,i),NSAMPLES(i),1e-8);
    %epsilon_rho_nnconvex(i) =getConfidence_nonconvex(SupportSize(i),1e-8,NSAMPLES(i));
    
    % Robustness
    beta=10^-8;
    OutWeJ = getWaitandJudgeEpsilon(SupportSize(1,i) ,NSAMPLES(i) ,beta);% get Wait-and-judge reliability
    epsilon_wait_and_judge =OutWeJ(end); % get Wait-and-judge reliability for the computed support size
    [EpsilonLU(1,i), EpsilonLU(2,i)] = epsLU(SupportSize(1,i),NSAMPLES(i), beta);
    Epsi =Get_Epsilon_Individual_Requirements(g_fun,delta_N,Results.dopt,beta,VaR); 
    Tempeps=Epsi';
    Epsilon_ind(:,i)=Tempeps(:);
 
end
% summarize resuts in a table

Designs=[dopt_exp2];
J=[Jopt];
CVaRpf=[CVAR_Fzero];
CVaR95=[CVAR95];
Pf_all=[Pf_all];
Pf_ind=[Pf_ind];
W_max=[Wmax];
G_max=[Gmax];
SN=[SupportSize];
%Epsilon=[epsilon_rho_convex];
Epsilon=[EpsilonLU];
%Table_2=[Designs;J;CVaRpf;CVaR95;Pf_all;Pf_ind;W_max;G_max;SN;Epsilon]
%%
Plot_Bounds_and_groundtruth(EpsilonLU,NSAMPLES,Pf_all)