%% MAIN
%% Mathematical setup
clc
clear variables
close all
% d \in \Theta a decision vector (eg. a system desing)
% \Theta= [LBd,UBd] a convex desing space
% delta = vector uncertain factors (scenario)
% DGM = the (unknown) data generating mechanism producing delta
% n_g= number of reliability requirements
% n_d = number of decision variables
% n_\delta = number of random variables
% J(d) = a convex cost function for taking a decision d
% g_j(d,\delta)>0 implies that the desing d fails to meet requirement j for the uncertain scenario \delta
% w(d,\delta)=\max\limits_{j\in\{1,..,n_g \}| g_j(d,\delta) %worst-case reliability performance function

%% CONSIDER THE FOLLOWING SCENARIO RBDO PRORAMS WITH SOFT-CONSTRAINTS:
%% Program-1
% min_{d\in \Theta , \zeta^{(i)}>0} \ lbrace J(d) +\rho \sum\limits_{i=1}^{N}  \zeta^{(i)}
% Such that: w(d,\delta{(i)} \leq \zeta^{(i)} \rbrace
% where $\rho$ is a parameter weighting the cost of violation for w and
% \zeta are slack variables used to relax reliability constraints for all requirements

%% Program-2
% min_{d\in \Theta , \zeta_j^{(i)}>0} \ lbrace J(d) +\sum\limits_{j=1}^{n_g} \rho_j \sum\limits_{i=1}^{N}  \zeta_j^{(i)}
% Such that: g_j(d,\delta) \leq \zeta_j^{(i)} i=1,...,N,~j=1,..,n_g\rbrace
% where $\rho_j$ are parameters weighting the cost of violation on the reliability requirement g_j
% \zeta_j are slack variables used to relax reliability constraints for requirement j

%% Robustness
% V(d*) = violation probability  i.e.  V(d*)=P[delta: w(d*,\delta)>0]
% epsilon = reliability parameter
% epsilon-robust solution V(d*)<epsilon
% beta = reliability parameter
% sN = number of support scenario


%% Load RBDO problem
CASESTUDY=1; % select case study
N=200;  % available samples
[g_fun,delta,dn,LBd,UBd,DGM,Nd,Ng]=Select_Convex_Case_Study(N,CASESTUDY);
N_testing=10^6;  % number of scenarios for validation
delta_testing=DGM(N_testing); % generate more scenarios for validation 
J_fun= @(d) sum(d+d.^2); % define cost function
%J_fun= @(d) -sum(d); % define cost function
% g_fun=@(d,delta) g_fun(d,delta)./(1+abs(g_fun(d,delta)))
w_fun=@(d,delta) max(g_fun(d,delta),[],2); % Worst-Case-Performance function

%% 
Rel_baseline=ComputeReliabilityPerformance(dn ,delta,g_fun);  % get reliability
J_baseline=J_fun(dn);
% Rel_dpreviouspaper=ComputeReliabilityPerformance([0.895,-0.2,0.14],delta,g_fun);
% J_dpreviouspaper=J_fun([0.895,-0.2,0.14]);
%% 1.1) Chance-Constrained Conditional Value at Risk (CVAR) problem
% min {J(d): CVaR_w(alpha) <= 0 }
% where CVaR_w(alpha) <= 0 ...this CVaR condition implies P[w>0] <= alpha
alpha=0.85; % define an alpha level
[Resultscvar,X_opt_cvar,J_opt_CVAR,~,~]=CVaR_CCP(delta,alpha,w_fun,J_fun,dn,LBd,UBd);
Rel_CVaR_opt=ComputeReliabilityPerformance(Resultscvar.dopt ,delta,g_fun);  % get reliability

%% R.T. Rockafellar 2010 RESS
NsamplesGMM=10^3;
Nd=length(dn); % get number of desing variables
LB=[LBd,-inf,zeros(1,NsamplesGMM)]; % lower bounds on d 
UB=[UBd,+inf,+inf*ones(1,NsamplesGMM)]; %upper bounds on d  
X0=[dn,10^4, zeros(1,NsamplesGMM)]; % an initial guess 
obj= @(x) J_fun(x(1:Nd)); % The cost function  
 
GMModel = fitgmdist(delta,size(delta,2)+3); % gaussian Mixture model 
 
% optionsGA = optimoptions('ga','PopulationSize',10000); % spq is needed to guaratnee high accuracy
%optionsGA
 % Scenario_Pf_constraint compute the constraints defined as w(d,delta_i)>\zeta_i 
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',50000,'Display','iter'); % spq is needed to guaratnee high accuracy
[X_opt,J_opt,exitflag,output]= fmincon(@(x) obj(x),... % objective function
             X0,[],[],[],[], LB, UB,...% initial guess and design bounds
              @(x) CVaR_Constraint_Rockafellar2010(x,GMModel,w_fun,alpha,Nd,NsamplesGMM),options); % non linear constraints and options
 Rel_CVaR_opt_rockafellar=ComputeReliabilityPerformance(X_opt(1:Nd) ,delta,g_fun);  % get reliability
 
%% 1.2) run scenario Program-1
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',50000,'Display','iter'); % spq is needed to guaratnee high accuracy
RHO= 100; % weight given to the violation
VaR=0; % this defines a value-at-risk for w (VaR=0 corresponds to a classical reliability problem)
[Results,X_opt,Obj_opt,~,~]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dn,LBd,UBd,options);
Rel_sp_opt=ComputeReliabilityPerformance(Results.dopt ,delta,g_fun);  % get reliability
sN=Results.Support.Size; % get size of the support sceanarios
sNj=sum(Rel_sp_opt.g>=0);
% 1.3) Robustness assessment
BetaConf=10^-8; % confidence level
OptType='convex_soft';
epsilon_soft_all=Epsilon_bounds(sN,N,BetaConf,OptType,Nd);% get Wait-and-judge epsilon
epsilon_wj=Epsilon_bounds(sN,N,BetaConf,'convex_w&j',Nd);% get Wait-and-judge epsilon
epsilon_nnconvex=Epsilon_bounds(sN,N,BetaConf,'nnconvex',Nd);% get Wait-and-judge epsilon
epsilon_soft_gj=Get_Epsilon_Individual_Requirements(g_fun,delta,Results.dopt,BetaConf,VaR);
Rel_sp_opt_trueprobability=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);  % get reliability

%% 2) Run Program-1 for increasing number of samples N=linspace(,,)

Nscenarios=[50,500,1000,1500,2000];nN_max=Nscenarios(end);

%Nscenarios=[50,500,1000,1500,2000]/10;nN_max=Nscenarios(end);
options.Display='off';
% prealocate memory
[CVAR95,CVAR_Fzero,Pf_all,Wmax,Jopt,SN,epsilon_rho_nnconvex,epsilon_wj]=deal(zeros(1,length(Nscenarios)));
dopt=zeros(Nd,length(Nscenarios));
[Gmax,Pf_ind]=deal(zeros(Ng,length(Nscenarios)));
EpsilonLU=zeros(2,length(Nscenarios));
Epsilon_ind=zeros(Ng*2,length(Nscenarios));

for i=1:length(Nscenarios) % loop for Nsamples
    rng default
    delta_N=DGM(Nscenarios(i));
    disp(num2str(Nscenarios(i)))
    [Results,X_opt,Obj_opt,exitflag,output]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta_N,J_fun,RHO,dn,LBd,UBd,options);
    dopt(:,i)=Results.dopt;
    Jopt(i)=Results.Jopt;
    SN(i)=Results.Support.Size;
    
    % Robustness analysis
    OptType='convex_soft';
    EpsilonLU(:,i)=Epsilon_bounds(SN(i),Nscenarios(i),BetaConf,OptType,Nd);% get Wait-and-judge epsilon
    epsilon_wj(i)=Epsilon_bounds(SN(i),Nscenarios(i),BetaConf,'convex_w&j',Nd);% get Wait-and-judge epsilon
    epsilon_rho_nnconvex(i)=Epsilon_bounds(SN(i),Nscenarios(i),BetaConf,'nnconvex',Nd);% get Wait-and-judge epsilon
    Temp=Get_Epsilon_Individual_Requirements(g_fun,delta_N,Results.dopt,BetaConf,VaR)';
    Epsilon_ind(:,i)=Temp(:);
    
    % get reliability  scores for the test set of scenarios
    Rel=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);
    CVAR_Fzero(i)=CVaR(Rel.w,1-Rel.Pf_all); CVAR95(i)=CVaR(Rel.w,0.95);
    Pf_all(i)=Rel.Pf_all; Pf_ind(:,i)=Rel.Pf_individual;
    Wmax(i)=Rel.w_max; Gmax(:,i)=Rel.g_max;
end
Plot_Bounds_and_groundtruth(EpsilonLU,Nscenarios,Pf_all)

%% 3) Run Program-1 for increasing number of samples N=linspace(,,)
%% and different Values-at-risk LambdaVal=linspace(,,)
Nlambda=6;
LambdaVal= linspace(-1.5,1.5,Nlambda);
% prealocate memory
[AlphaValueAtRisk,CVAR_Fzero,Pf_all,Jopt,SN_lambda,...
    epsilon_rho_convex,epsilon_lambda_nnconvex,epsilon_lambda_wj,...
    Alpha_true,CVAR_true,Pf_true,epsilon_true_all]=deal(zeros(1,Nlambda));
dopt=zeros(Nlambda,Nd);
[EpsilonLU_lambda]=deal(zeros(2,Nlambda));
epsilon_true_ind=zeros(Ng,Nlambda);
[Epsilon_ind_lambda]=deal(zeros(Ng*2,Nlambda));
for nN= sort(Nscenarios,'descend') % loop different N values
    rng default
    delta_N=DGM(nN); % sample nN scenarios
    for i=1:Nlambda % loop different lambda (values at risk)
        % lambda = 0 classical RBOD (g>0 is violation)
        % lambda > 0  less stringent requirement on g (g-VaR>0 is violation)
        % lambda < 0  harden the requirements on g (g+VaR>0 is violation)
        disp(['Optimizing for: VaR =' num2str(LambdaVal(i)) ' and N =' num2str(nN)]);
        VaR=LambdaVal(i); % this defines the value at risk to compute the CVAR
        [Results,X_opt,Obj_opt,exitflag,output]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta_N,J_fun,RHO,dn,LBd,UBd,options);
        AlphaValueAtRisk(i)=Results.AlphaValueAtRisk;
        dopt(i,:)=Results.dopt;
        Jopt(i)=Results.Jopt;
        SN_lambda(i)=Results.Support.Size;
        
        % Robustness analysis
        OptType='convex_soft';
        EpsilonLU_lambda(:,i)=Epsilon_bounds(SN_lambda(i),nN,BetaConf,OptType,Nd);% get Wait-and-judge epsilon
        epsilon_lambda_wj(i)=Epsilon_bounds(SN_lambda(i),nN,BetaConf,'convex_w&j',Nd);% get Wait-and-judge epsilon
        epsilon_lambda_nnconvex(i)=Epsilon_bounds(SN_lambda(i),nN,BetaConf,'nnconvex',Nd);% get Wait-and-judge epsilon
        Temp=Get_Epsilon_Individual_Requirements(g_fun,delta_N,Results.dopt,BetaConf,VaR)';
        Epsilon_ind_lambda(:,i)=Temp(:);
        
        %% compute reliability for the optimized desing a
        Rel=ComputeReliabilityPerformance(Results.dopt,delta_N,g_fun);  % get reliability
        CVaR_lambda(i)=mean(Rel.w(Rel.w>=VaR));
        CVAR_Fzero(i)=CVaR(Rel.w,1-Rel.Pf_all);
        CVAR95(i)=CVaR(Rel.w,0.95);
        W_max(i)=Rel.w_max;
        Percentile_w(i,:)=prctile(Rel.w,[85:1:100]);
        Pf_all(i)=Rel.Pf_all;
        
        %% true statistics on the larger dataset (10^5 scenarios)
        res1e5=ComputeReliabilityPerformance(Results.dopt,delta_testing,g_fun);
        epsilon_true_all(i)=mean(res1e5.w>=VaR); % true violation probability
        CVAR_true(i)=mean(res1e5.w(res1e5.w>=VaR));
        CVAR_95true(i)=CVaR(res1e5.w,0.95);
        W_max_true(i)=res1e5.w_max;
        CVAR_Fzerotrue(i)=CVaR(res1e5.w,1-res1e5.Pf_all);
        Pf_true(i)=mean(res1e5.w>=0);
        epsilon_true_ind(:,i)=mean(res1e5.g>=VaR);
    end
    %% plot joint requirements bounds
    figure(1)
    hold on;
    if nN == Nscenarios(1)
        patch([LambdaVal fliplr(LambdaVal)], [EpsilonLU_lambda(1,1:Nlambda) fliplr(EpsilonLU_lambda(2,1:Nlambda))], 'b','Facealpha',0.08)
        %  plot(LambdaVal,epsilon_true_all,'k','LineWidth',2)
    elseif  nN == Nscenarios(end)
        patch([LambdaVal fliplr(LambdaVal)], [EpsilonLU_lambda(1,1:Nlambda) fliplr(EpsilonLU_lambda(2,1:Nlambda))], 'b','Facealpha',0.08)
        plot(LambdaVal,epsilon_true_all,'-or','LineWidth',2)
     %   legend('$[\underline{\epsilon}(s_N^\star),\overline{\epsilon}(s_N^\star)]$','${P}[w(d,\delta)>\lambda]$','Interpreter','latex')
    else
        patch([LambdaVal fliplr(LambdaVal)], [EpsilonLU_lambda(1,1:Nlambda) fliplr(EpsilonLU_lambda(2,1:Nlambda))], 'b','Facealpha',0.08)
    end
    grid on; box on;
    xlim([min(LambdaVal),max(LambdaVal)])
    ylim([0,1])
    xlabel('\lambda');
    ylabel('$P[w>\lambda]$','Interpreter','latex')
    set(gca,'FontSize',18)
    %% Plot individual requirements bounds
    figure(2)
    for i=1:Ng
        EPSBOUND=Epsilon_ind_lambda([2*(i-1)+1,2*(i-1)+2],:);
        subplot(ceil(Ng/2),2,i);
        hold on
        if nN ==  Nscenarios(1)
            patch([LambdaVal fliplr(LambdaVal)], [EPSBOUND(1,1:Nlambda) fliplr(EPSBOUND(2,1:Nlambda))], 'b','Facealpha',0.1)
            %   plot(LambdaVal,epsilon_true_ind(i,:),'k','LineWidth',2)
        elseif  nN == Nscenarios(end)
            patch([LambdaVal fliplr(LambdaVal)], [EPSBOUND(1,1:Nlambda) fliplr(EPSBOUND(2,1:Nlambda))], 'b','Facealpha',0.1)
            plot(LambdaVal,epsilon_true_ind(i,:),'-or','LineWidth',2)
        else
            patch([LambdaVal fliplr(LambdaVal)], [EPSBOUND(1,1:Nlambda) fliplr(EPSBOUND(2,1:Nlambda))], 'b','Facealpha',0.1)
        end
        grid on; box on;
        xlim([min(LambdaVal),max(LambdaVal)])
        ylim([0,1])
        xlabel('\lambda')
        ylabel(['$P[g_' num2str(i) '>\lambda]$'],'Interpreter','latex')  
        set(gca,'FontSize',18) 
    end
    
end
figure(1)
        legend('$[\underline{\epsilon}(s_N^\star),\overline{\epsilon}(s_N^\star)]$','${P}[w(d^\star,\delta)\geq\lambda]$ ($\mathcal{D}_{10^6}$)','Interpreter','latex')
figure(2)
    for i=1:Ng 
        subplot(ceil(Ng/2),2,i);
        legend(['$[\underline{\epsilon}(s_N^\star),\overline{\epsilon}(s_N^\star)]$'],['${P}[g_' num2str(i) '(d^\star,\delta)\geq\lambda]$ ($\mathcal{D}_{10^6}$)'],'Interpreter','latex')
    end
%% 3) Run Program-1 for incresing cost of violation Rho = linspace(,,)
Nrho=50
VaR=0;% value at risk (violation is for g>=VaR)
RHOLINSPACE=linspace(0.001,1,Nrho);
%RHOLINSPACE=[0.001,0.01,0.05,0.1,0.5,1];
%prealocate memory
[SN,Cost,Pf_all_rho]=deal(zeros(1,Nrho));
d_opt=zeros(Nd,Nrho);
EpsilonLU_rho=deal(zeros(2,Nrho));
[Pf_ind_rho]=deal(zeros(Ng,Nrho));
epsilon_rho_wj=zeros(1,Nrho);
Epsilon_ind_rho=zeros(2*Ng,Nrho);
for i=1:Nrho
    % weight given to the violation
    RHO= RHOLINSPACE(i);
    disp(['Optimizing for: rho =' num2str(RHOLINSPACE(i)) ]);
    
    % run optimizer
    [Results,X_opt,J_opt,exitflag,~]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dn,LBd,UBd,options);
    d_opt(:,i)=Results.dopt;
    SN(i)=Results.Support.Size; % get size of the support sceanarios
    Cost(i)=Results.Jopt;
    
    % Robustness analysis
 
   % epsilon_rho_wj(i)=Epsilon_bounds(SN(i),N,BetaConf,'convex_w&j',Nd);% get Wait-and-judge epsilon
   % epsilon_rho_nnconvex(i)=Epsilon_bounds(SN(i),N,BetaConf,'nnconvex',Nd);% get Wait-and-judge epsilon
    Temp=Get_Epsilon_Individual_Requirements(g_fun,delta,Results.dopt,BetaConf,VaR)';
    Epsilon_ind_rho(:,i)=Temp(:);
    
    % get reliability performance for the delta_testing
    Rel_sp_opt=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);
    Pf_all_rho(i)=Rel_sp_opt.Pf_all;
    Pf_ind_rho(:,i)=Rel_sp_opt.Pf_individual;
    %wmax_rho(i)=Rel_sp_opt.w_max;
end

Plot_Bounds_and_groundtruth(EpsilonLU_rho,RHOLINSPACE,Pf_all_rho)
%% 4) Run Scenario Program-2 for MeshGrid of (rho_1,rho_2) check trade off between individual Pf
%  Program-2
% min_{d\in \Theta , \zeta_j^{(i)}>0} \ lbrace J(d) +\sum\limits_{j=1}^{n_g} \rho_j \sum\limits_{i=1}^{N}  \zeta_j^{(i)}
% Such that: g_j(d,\delta) \leq \zeta_j^{(i)} i=1,...,N,~j=1,..,n_g\rbrace
% where $\rho_j$ are parameters weighting the cost of violation on the reliability requirement g_j

if Ng<=1; error('Program 2 should be used when there are 2 or more reliability performance function gj'); end

GIdxi=1; % i th requirement index to be analysed
Gidxj=2; % j th requirement index to be analysed
Ngridpoints=50; % number of gird points in the i,j space
RHOivec=linspace(10^-3, 0.1 , Ngridpoints);
RHOjvec=linspace(10^-3, 0.1 , Ngridpoints);
[RHOi,RHOj]=meshgrid(RHOivec,RHOjvec);
RHOs= [RHOi(:) RHOj(:)]; % vector of rho_ij costs of vioaltion coordinates
Nrhos=size(RHOs,1);
VaR=0; % value at risk (violation is for g>=VaR)
%prealocate memory
[SNj,Pf_ind]=deal(zeros(Nrhos,Ng));
Epsilon_ind_multirho=zeros(2*Ng,Nrhos);
[Pf_all,Jopt]=deal(zeros(1,Nrhos));
dopt=zeros(Nd,Nrhos);

Rhoj=ones(1,Ng); % cost vector
for rj=1:Nrhos
    
    % get cost of violation
    Rhoj([GIdxi,Gidxj])=RHOs(rj,:);
    disp(['Solving optimization ' num2str(rj) '/' num2str(Nrhos) '  for costs: ' num2str(Rhoj) ])
    
    % run optimizer
%     g_norm =  @(d,delta) g_fun(d,delta)./(1+abs(g_fun(d,delta)));
%     [Results,X_opt,J_opt,exitflag,output]=...
%         ScenarioOptimizerCVaR_program2(VaR,g_norm,delta,J_fun,Rhoj,dn,LBd,UBd,options);
    [Results,X_opt,J_opt,exitflag,output]=...
       ScenarioOptimizerCVaR_program2(VaR,g_fun,delta,J_fun,Rhoj,dn,LBd,UBd,options);
    
    SNj(rj,:)=Results.Support.Size;
    Jopt(rj)=J_fun(X_opt(1:Nd));
    dopt(:,rj)= X_opt(1:Nd);
    
    % Robustness analysis
    OptType='convex_soft';
    for i=1:length(SNj(rj,:))
        %EpsilonLU_lambda(:,rj)=Epsilon_bounds(SN_lambda(i),nN,BetaConf,OptType,Nd);% get Wait-and-judge epsilon
        epsilon_lambda_wj(rj,i)=Epsilon_bounds(SNj(rj,i),N,BetaConf,'convex_w&j',Nd);% get Wait-and-judge epsilon
        epsilon_lambda_nnconvex(rj,i)=Epsilon_bounds(SNj(rj,i),N,BetaConf,'nnconvex',Nd);% get Wait-and-judge epsilon
    end
    Temp =Get_Epsilon_Individual_Requirements(g_fun,delta,X_opt(1:Nd),BetaConf,VaR)';
    Epsilon_ind_multirho(:,rj)=Temp(:);
    
    %  Reliability verification
    Rel_multiRHOG=ComputeReliabilityPerformance(Results.dopt,delta,g_fun);
    Pf_all(rj)=Rel_multiRHOG.Pf_all;
    Pf_ind(rj,:)=Rel_multiRHOG.Pf_individual;
end

figure
ThresholdPf1=0.2;
ThresholdPf2=0.5;
ThresholdCost=10;
IDX_toplot= Jopt<ThresholdCost & Pf_ind(:,1)'<ThresholdPf1 & Pf_ind(:,2)'<ThresholdPf2;
surf(RHOi,RHOj,reshape(Jopt.*IDX_toplot,Ngridpoints,Ngridpoints),'FaceColor','black','FaceAlpha',0.3,'EdgeColor',[0.6 0.6 0.6],'EdgeAlpha',0.3)
xlabel('$\rho_i$','Interpreter','latex')
ylabel('$\rho_j$','Interpreter','latex')
set(gca,'FontSize',20)
hold on
TEMP=Pf_ind(:,1).*IDX_toplot';
TEMP(TEMP==0)=inf;
scatter3(RHOs(min(TEMP)==TEMP,1),RHOs(min(TEMP)==TEMP,2),20,'sb','filled') 
TEMP=Pf_ind(:,2).*IDX_toplot';
TEMP(TEMP==0)=inf;
scatter3(RHOs(min(TEMP)==TEMP,1),RHOs(min(TEMP)==TEMP,2),20,'ob','filled')
TEMP=Pf_all.*IDX_toplot;
TEMP(TEMP==0)=inf;
scatter3(RHOs(min(TEMP)==TEMP,1),RHOs(min(TEMP)==TEMP,2),20,'hk','filled')
TEMP=Jopt.*IDX_toplot;
TEMP(TEMP==0)=inf;
scatter3(RHOs(min(TEMP)==TEMP,1),RHOs(min(TEMP)==TEMP,2),20,'dk','filled')


figure
scatter3(Pf_ind(:,1),Pf_ind(:,2), Jopt')
[ p, idxs] = paretoFront( -[Pf_ind Jopt'] );
hold on; scatter3(Pf_ind(idxs,1),Pf_ind(idxs,2), Jopt(idxs)','+r')
figure;
IDX=zeros([size(RHOs,1),1]);
IDX(idxs)=1;
gplotmatrix([RHOs SNj Jopt'],[],IDX,'br')
% plot
figure
surf(RHOi,RHOj,reshape(SNj(:,1)./(SNj(:,1)+SNj(:,2)),Ngridpoints,Ngridpoints),'FaceColor','blue','FaceAlpha',0.3,'EdgeColor',[0.1 0.1 0.9],'EdgeAlpha',0.3)

figure
subplot(1,2,1)
surf(RHOi,RHOj,reshape(Pf_ind(:,1),Ngridpoints,Ngridpoints),'FaceColor','black','FaceAlpha',0.3,'EdgeColor',[0.6 0.6 0.6],'EdgeAlpha',0.3)
hold on
surf(RHOi,RHOj,reshape(Epsilon_ind_multirho(1,:),Ngridpoints,Ngridpoints),'FaceColor','red','FaceAlpha',0.1,'EdgeColor',[0.9 0 0],'EdgeAlpha',0.1)
surf(RHOi,RHOj,reshape(Epsilon_ind_multirho(2,:),Ngridpoints,Ngridpoints),'FaceColor','red','FaceAlpha',0.1,'EdgeColor',[0.9 0 0],'EdgeAlpha',0.1)
surf(RHOi,RHOj,reshape(Normalizein01(Jopt),Ngridpoints,Ngridpoints),'FaceColor','blue','FaceAlpha',0.3,'EdgeColor',[0.1 0.1 0.9],'EdgeAlpha',0.3)
xlabel('$\rho_i$','Interpreter','latex')
ylabel('$\rho_j$','Interpreter','latex')
zlabel('$\frac{s_{N,i}}{s_{N,i}+s_{N,j}}$','Interpreter','latex')

subplot(1,2,2)
surf(RHOi,RHOj,reshape(Pf_ind(:,2),Ngridpoints,Ngridpoints),'FaceColor','black','FaceAlpha',0.3,'EdgeColor',[0.6 0.6 0.6],'EdgeAlpha',0.3)
hold on
surf(RHOi,RHOj,reshape(Epsilon_ind_multirho(3,:),Ngridpoints,Ngridpoints),'FaceColor','red','FaceAlpha',0.1,'EdgeColor',[0.9 0 0],'EdgeAlpha',0.2)
surf(RHOi,RHOj,reshape(Epsilon_ind_multirho(4,:),Ngridpoints,Ngridpoints),'FaceColor','red','FaceAlpha',0.1,'EdgeColor',[0.9 0 0],'EdgeAlpha',0.2)
surf(RHOi,RHOj,reshape(Normalizein01(Jopt),Ngridpoints,Ngridpoints),'FaceColor','blue','FaceAlpha',0.3,'EdgeColor',[0.1 0.1 0.9],'EdgeAlpha',0.3)

xlabel('$\rho_i$','Interpreter','latex')
ylabel('$\rho_j$','Interpreter','latex')
zlabel('$\frac{s_{N,i}}{s_{N,i}+s_{N,j}}$','Interpreter','latex')
%set(gca,'FontSize',22)
% for i=1:Ng
%     subplot(ceil(Ng/2),2,i);
%     Plot_Bounds_and_groundtruth(Epsilon_ind_multirho([2*(i-1)+1,2*(i-1)+2],:),1:1:10,Pf_ind(:,i))
% end
