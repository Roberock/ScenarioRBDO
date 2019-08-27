function [ScenarioRobustness]=Support_Constraints_addMethod_SP3(obj,Theta_opt,alpha,Theta0,Tol,G2ex)

%% FindSupportConstraints
% This function extract support contraints for a solution d of a
% non-convex scenario program
%  THE PROCEDURE IS SIMILAR TO THE ONE IN: Campi et al. 2018, A General Scenario Theory for Nonconvex Optimization and Decision Making
%% THE PROCEDURE WORKS AS FOLLOWS:
% 1) Set the initial solution to the original inital guess $\theta_{0}=(d_0,\alpha_0)$;
% 2) Initialize an empty set of support constraint $\mathcal{S}\leftarrow \emptyset$; (IncludedIndices)
% 3) For each $\delta\in \mathcal{D_\delta}$, $j=1,..,n_g$ and  in correspondence of $\theta_{old}$ compute $g_j(d,\delta)$;
% 4) Find worst scenario constrain $i=\arg\max\limits_{i,j} ~~g_j(d,\delta^{(i)})-\alpha_j$ and update the support $\mathcal{S}\leftarrow \mathcal{S} \cup \lbrace i \rbrace$ and find optimal solution $\theta_{new}=SP(\mathcal{S})$ with starting solution $(d_0,\alpha_0)$
% 5) If $\theta_{new}=(d^*,\alpha^*)$ stop procedure and return the support set cardinality $s_N^*=|\mathcal{S}|$;. Otherwise, set $\theta_{old}=\theta_{new}$ and go to step (3);
%% INPUT
% d_opt=optimal design for the scenario program
% delta= matrix of scenarios [N x ndelta], observations from an unknown probabilistic model P
% ScenarioOptData
% e.g. ScenarioOptData.alpha= objective function handle for the scenario optimization
% alpha=@(x,Nalpha,Nd) sum((x(Nd+1:Nd+Nalpha)>=0).*x(Nd+1:Nd+Nalpha)); % objective function
% SysData = sytem solver data
% e.g. SysData.MyminfunG= performance function handle
% Nalpha= objective function number of alpha
%% OUTPUT
% Cardinality= cardinality of the support set;
% IncludedIndices = indices of the scenarios which are in the support set
% optionsFMINCON=optimoptions('fmincon','Display','iter');
% Error= relative error of the newly computed solution w.r.t. the original problem
%% Scenario optimization and system solver data
delta=obj.Sub_delta;  % The scenarios examined
fun_g=obj.g_fun;   % performance function handle
NSubdelta=obj.NSubdelta; % number of scenarios examined
Nrv=obj.Nrv; % scenarios dimensionality
%% 2) 2) Initialize an empty set of support constraint and prealocate memory
Gdeltaidx=[];
RMSE=zeros(1,NSubdelta); 
CountIter=1; % count the iterations 
D_now=Theta_opt;
delta_iter=[delta]; % add a column with scenario indices
Gdeltaidx=[];
RemainingIndx=1:1:NSubdelta;
%% iterativelly add Scenarios 
tic
while true
    % 1) for each scenario (remove index) and for desing Dnow, compute the reliability performance G
    G= fun_g(delta_iter(:,1:Nrv),D_now(1:obj.Nd));
    %% Identify the the worst case scenario WCS and update L
    % normalize g values and compute W
    gnorm=G./(1+abs(G));%% normalize G values in [-1 1] 
    LogicIdx= find(gnorm(:,G2ex)==max(gnorm(RemainingIndx,G2ex))); % scenario index, pick scenarios delta: W(delta)=max(W) or
    % LogicIdx=W>=prctile(W,98);  % scenario index, pick scenarios with W>98 percentile 
    RemainingIndx=setdiff(RemainingIndx,LogicIdx);
    Gdeltaidx=union(Gdeltaidx,LogicIdx); % include new worst-case in the worst case indices
   
    %% 2) run the optimization with L as constraint set
   % obj.Sub_delta=L;  
    % [D_now,Gamma_min,exitflag,output]= obj.Scenario_RBDO(alpha, Theta0);
    [D_now,~,~,~]= obj.Optimize_SP3(alpha, Theta0,G2ex,Gdeltaidx); 
   % Error= sqrt((D_now- Theta_opt).^2)./Theta_opt; % Check Error 
   % RMSE(CountIter)= sqrt(mean((D_now- Theta_opt).^2)); 
     RMSE(CountIter)= sqrt(mean((D_now([1:obj.Nd,G2ex])- Theta_opt([1:obj.Nd,G2ex])).^2)); 
    % if the solution is equal to the original one with a small error, break the while
    if   RMSE(CountIter)<Tol || CountIter>=NSubdelta %GammaError<1e-3 && RMSE<1e-3
        break
    end 
    CountIter=CountIter+1; % count iter number
end

%% Now reduce to a set with min cardinality
if length(Gdeltaidx)==1
    
    ScenarioRobustness.Support_Set_Min=Gdeltaidx;
    ScenarioRobustness.Support_Set_nonMin=Gdeltaidx;
    ScenarioRobustness.Cardinality=1;
    ScenarioRobustness.CountIter=CountIter+0;
    ScenarioRobustness.ComputationalTime.ConstranScenarioAddMethod=toc;
    ScenarioRobustness.RMSE=RMSE(1:CountIter);
    
else

    ScenarioRobustness.ComputationalTime.ConstranScenarioAddMethod=toc;
    ResultsPart2=obj.ScenarioConstraintsSp3_removeMethod(G2ex,Tol,alpha,Theta0,Gdeltaidx);
    ScenarioRobustness.ComputationalTime.ConstranScenarioRemoveMethod=toc;
    ScenarioRobustness.Support_Set_Min=ResultsPart2.Support_Set;
    ScenarioRobustness.Support_Set_nonMin=Gdeltaidx;
    ScenarioRobustness.Cardinality=ResultsPart2.Cardinality;
    ScenarioRobustness.CountIter=CountIter+ResultsPart2.CountIter;
    %ScenarioRobustness.ComputationalTime=toc;
    ScenarioRobustness.RMSE=RMSE(1:CountIter);
end
end