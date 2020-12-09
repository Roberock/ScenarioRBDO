function [Results,X_opt,J_opt,exitflag,output]=ScenarioOptimizerCVaR_program2(VaR,g_fun,delta,J_fun,Rho,dn,LBd,UBd,options)
 %% Input
 % VaR value-at-risk (lamba) for each g_j
 % g_fun = @(d,delta) .... the reliability performance function vectorized
 % delta = matrix containing the scenarios [N x Ndelta]
 % J_fun = cost function for the decision d
 % Rho = vector of weights...the costs of reliability violations for g_j
 % dn = nominal desing (initial guess in the optimizer)
 % LBd,UBd, lower and upper bounds for the desing 
 % options optimization options
 
%% define the optimization program
%   min_{d,\zeta_i}       Obj(d)
%          s.t.           Obj(d)=J(d) +\sum\limits_{j=1}^{ng}\rho_j\sum\limits_{i=1}^{N} \zeta_{j,i}
%          s.t.           g_j(d,\delta{(i)} <= \zeta_{j,i},~ j=1,...,ng, i=1,...,N
%          s.t.           d \in \Theta & \zeta_{j,i} > 0
% where J(d) is a cost function
% $\zeta_{j,i}$ is the slack variable associated to the scenario i and requirement j 
% \sum\limits_{i=1}^{N} \zeta_{j,i} with \zeta_{j,i}>0 is proportional to the mean of g_j in the failure region
% $\rho_j$ is a scaling parameter (cost of violation for g_j)
% g_j(d,\delta) is the reliability performance fun. wher eg>0 implies failure
 
%% Start
Nd=length(dn); % get number of desing variables
N=size(delta,1);
Ng=length(g_fun(dn,delta(1,:))); 
LB=LBd; UB=UBd; % buld bounds for the optimization problem
for i=1:Ng
    LB=[LB, VaR*ones(1,N)]; % lower bounds on d and the slack variables \zeta_ij
    UB=[UB, inf*ones(1,N)]; %upper bounds on d and the slack variables \zeta_ij
end
X0=[dn, LB(Nd+1:end)]; % initial guess
% The cost function with slack variables J(d)+\sum_j rho_j sum_i \zeta_ij 
obj= @(x) J_fun(x(1:Nd)); % build objective function
for i=1:Ng
        obj=  @(x) obj(x)+  Rho(i)*sum(x(Nd+1+(i-1)*N:Nd+N+(i-1)*N)+LB(Nd+1+(i-1)*N:Nd+N+(i-1)*N));
end
 %% run fmincon optimzer 
 options.Algorithm='interior-point'; % generally we have a large scale problem 
 options.Display='off';
 % g_norm =  @(d,delta) g_fun(d,delta)./max(g_fun(dn,delta)); % CONSIDER THIS TO MAKE THE REQUIREMENTS COMPARABLE
 
 [X_opt,J_opt,exitflag,output]= fmincon(@(x) obj(x),... % objective function
        X0,[],[],[],[], LB, UB,...% initial guess and design bounds
        @(x) Scenario_Pf_constraint_multipleSoftG(x,delta,g_fun,Nd),options); % non linear constraints and options
    
 %% Collect Results
Results.Jopt = J_fun(X_opt(1:Nd));
Results.dopt = X_opt(1:Nd);
Results.Zopt = reshape(X_opt(Nd+1:end),N,Ng);  
Results.AlphaValueAtRisk = 1-(mean(Results.Zopt>VaR )+0.5/N); %alpha level corresponding to the VaR selected (add half step of the ECDF)
%% Support scenarios
Results.Support.Size = sum(g_fun(X_opt(1:Nd),delta)>=VaR);
Results.Support.Scenarios = (g_fun(X_opt(1:Nd),delta)>=VaR); 
end