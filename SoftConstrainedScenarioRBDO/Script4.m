Nlambda=10;
nN_max=1000;
[g_fun,delta,dn,LBd,UBd,DGM,Nd,Ng]=Select_Convex_Case_Study(nN_max,CASESTUDY);

RHO= 100; % weight given to the violation (Rho=0 means we have no gain from violating constraints)
LambdaVal= linspace(-2,2,Nlambda) ;
[AlphaValueAtRisk,CVAR_Fzero,Pf_all,Jopt,SupportSize,epsilon_rho_convex,epsilon_rho_nnconvex,Alpha_true,CVAR_true,Pf_true]=deal(zeros(1,Nlambda));
for nN= 100:100:nN_max
    deltaN=delta(1:nN,:);
    % for nN= [100, 1000, 2000]
    %deltaN=delta_testing(1:nN,:);
    
    for i=1:length(LambdaVal)
        
        display(num2str(i));
        VaR=LambdaVal(i); % this defines the value at risk to compute the CVAR
        [Results,X_opt,J_opt,exitflag,output]=ScenarioOptimizerCVaR(VaR,w_fun,deltaN,J_fun,RHO,dn,LBd,UBd,options);
        Rel=ComputeReliabilityPerformance(Results.dopt ,deltaN,g_fun);  % get reliability
        AlphaValueAtRisk(i)=Results.AlphaValueAtRisk;
        dopt_exp2(i,:)=Results.dopt;
        
        CVaR_lambda(i)=mean(Rel.w(Rel.w>=VaR));
        CVAR_Fzero(i)=CVaR(Rel.w,1-Rel.Pf_all);
        CVAR95(i)=CVaR(Rel.w,0.95);
        
        W_max(i)=Rel.w_max;
        Percentile_w(i,:)=prctile(Rel.w,[85:1:100]);
        
        Pf_all(i)=Rel.Pf_all;
        Jopt(i)=Results.Jopt;
        SupportSize(i)=Results.Support.Size;
        % examin scenario prescpective reliability
        beta=10^-8;
        OutWeJ =getWaitandJudgeEpsilon(SupportSize(i),nN,beta);
        epsilon_rho_convex(i)=OutWeJ(end);
        epsilon_rho_nnconvex(i) =getConfidence_nonconvex(SupportSize(i),beta,nN);
        [epsL, epsU] = epsLU(SupportSize(i),nN,beta);
        Epsilon_relaxation(:,i)=[epsL, epsU];
        %robustness
        
        Epsi =Get_Epsilon_Individual_Requirements(g_fun,deltaN,dopt_exp2(i,:),beta,VaR);
        Tempeps=Epsi';
        Epsilon_individual(:,i)=Tempeps(:);
        %% true statistics on the larger dataset (10^5 scenarios)
        
        res1e5=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);
        epsilon_true(i)=mean(res1e5.w>=VaR); % true violation probability
        CVAR_true(i)=mean(res1e5.w(res1e5.w>=VaR));
        CVAR_95true(i)=CVaR(res1e5.w,0.95);
        W_max_true(i)=res1e5.w_max;
        CVAR_Fzerotrue(i)=CVaR(res1e5.w,1-res1e5.Pf_all);
        Pf_true(i)=mean(res1e5.w>=0);
        epsilon_true_ind(:,i)=mean(res1e5.g>=VaR);
    end
    %% plot joint requirements bounds
    figure(1)
    patch([LambdaVal fliplr(LambdaVal)], [Epsilon_relaxation(1,1:length(LambdaVal)) fliplr(Epsilon_relaxation(2,1:length(LambdaVal)))], 'r','Facealpha',nN./(nN_max*10))
    grid on; box on; hold on; xlabel('\lambda');
    ylabel('$P[w>\lambda]$','Interpreter','latex')
    if nN == nN_max
        plot(LambdaVal,epsilon_true,'k')
        legend('$[\underline{\epsilon}(s_N^\star),\overline{\epsilon}(s_N^\star)]$','${P}[w(d,\delta)>\lambda]$','Interpreter','latex')
    end
    %% Plot individual requirements bounds
    figure(2)
    for i=1:Ng
        EPSBOUND=Epsilon_individual([2*(i-1)+1,2*(i-1)+2],:);
        subplot(ceil(Ng/2),2,i);
        patch([LambdaVal fliplr(LambdaVal)], [EPSBOUND(1,1:length(LambdaVal)) fliplr(EPSBOUND(2,1:length(LambdaVal)))], 'r','Facealpha',0.1)
        grid on
        box on
        xlabel('\lambda')
        ylabel(['$P[g_' num2str(i) '>\lambda]$'],'Interpreter','latex')
        hold on
        if nN == nN_max 
        plot(LambdaVal,epsilon_true_ind(i,:),'k','LineWidth',2)
         legend(['[\underline{\epsilon}(s_N^\star),\overline{\epsilon}(s_N^\star)]$','${P}[g_' num2str(i) '(d,\delta)>\lambda]$'],'Interpreter','latex')
         set(gca,'FontSize',18)
         end
        xlim([min(LambdaVal),max(LambdaVal)])
    end
    
end
%Idx2Keep=[1,4,8,12,16,20];
Idx2Keep=[1:1:6];
Table_3=[LambdaVal(Idx2Keep)
    dopt_exp2(Idx2Keep,:)'
    Jopt(Idx2Keep)
    CVaR_lambda(Idx2Keep)
    %CVAR_Fzero(Idx2Keep)
    CVAR95(Idx2Keep)
    Pf_all(Idx2Keep)
    W_max(Idx2Keep)
    SupportSize(Idx2Keep)
    Epsilon_relaxation(:,Idx2Keep)  ]
Table_3_verification=[CVAR_true(Idx2Keep);CVAR_Fzerotrue(Idx2Keep);...
    CVAR_95true(Idx2Keep);Pf_true(Idx2Keep);W_max_true(Idx2Keep);epsilon_true(Idx2Keep)]