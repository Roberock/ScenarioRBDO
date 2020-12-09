function [Results,X_opt,J_opt,exitflag,output]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,Rho,dn,LBd,UBd,options)
 %% Input
 % VaR value-at-risk (lamba)  
 % w_fun = (d,delta) .... worst-case reliability performance function
 % delta = matrix containing the scenarios [N x Ndelta]
 % J_fun = cost function for the decision d
 % Rho = scalar weighting the cost of reliability violations
 % dn = nominal desing (initial guess in the optimizer)
 % LBd,UBd, lower and upper bounds for the desing 
 % options optimization options
 
%% define the optimization program
%   min_{d,\zeta_i}       Obj(d)
%          s.t.           Obj(d)=J(d) +\rho \sum\limits_{i=1}^{N} \zeta_i
%          s.t.           w(d,\delta{(i)} <= \zeta_i
%          s.t.           d\in \Theta & \zeta_i>VaR
% where J(d) is a cost function
% $\zeta_j^{(i)}$ are slack variables 
% \sum\limits_{i=1}^{N} \zeta_i with \zeta_i>VaR is proportional to CVAR
% $\rho_j$ is a scaling parameter (importance of CVAR vs J(d))
% and w(d,\delta)=max_j g_j(d,\delta) is the worst-case performance fun.
%% Example
% delta= [7.79,-9.59;-0.16,4.79;-4.25,-1.72;-1.13,-1.69]; % [Nsamples x Ndelta]
%  g1= @(d,delta)  -delta(:,1).*d(1) + d(2) +delta(:,2) ;   % linear Performance Functions 1
%  g2= @(d,delta)  -delta(:,2).*d(2) + d(1) +delta(:,1);  % linear Performance Functions 2
%  g_fun=@(d,delta)[g1(d,delta), g2(d,delta) ]; % Performance Functions
%  w_fun=@(d,delta) max(g_fun(d,delta),[],2); % Worst-Case-Performances Across Samples (normalized [-1 1]) 
% LBd=[-5,-5]; UBd=[5,5]; 
% dn=[1,1]
% J_fun= @(d) -sum(d); % first part of the cost function J(d)  
% RHO= 100; % weight given to the violation 
% VaR=0; % this defines the value at risk to compute the CVAR 
%% Start
Nd=length(dn); % get number of desing variables
Ndelta=size(delta,1);
LB=[LBd, VaR*ones(1,Ndelta)]; % lower bounds on d and the slack variables \zeta_i
UB=[UBd, inf*ones(1,Ndelta)]; %upper bounds on d and the slack variables \zeta_i
X0=[dn, LB(Nd+1:end)]; % an initial guess 
obj= @(x) J_fun(x(1:Nd))+ Rho*sum(x(Nd+1:end)+LB(Nd+1:end)); % The cost function with slack variables J(d)+sum_i \zeta_i 
 %% run fmincon optimzer 
 % Scenario_Pf_constraint compute the constraints defined as w(d,delta_i)>\zeta_i
[X_opt,J_opt,exitflag,output]= fmincon(@(x) obj(x),... % objective function
               X0,[],[],[],[], LB, UB,...% initial guess and design bounds
              @(x) Scenario_Pf_constraint(x,delta,w_fun,Nd),options); % non linear constraints and options
  
 %% Collect Results
Results.Jopt = J_fun(X_opt(1:Nd));
Results.dopt = X_opt(1:Nd);
Results.Zopt = X_opt(Nd+1:end);  
Results.AlphaValueAtRisk = 1-(mean(Results.Zopt>VaR )+0.5/Ndelta); %alpha level corresponding to the VaR selected (add half step of the ECDF)
%% Support scenarios
Results.Support.Size = sum(w_fun(X_opt(1:Nd),delta)>=VaR);
Results.Support.Scenarios = (w_fun(X_opt(1:Nd),delta)>=VaR); 
end