function [Results,X_opt,J_opt,exitflag,output]=CVaR_CCP(delta,alpha,w_fun,J_fun,dn,LBd,UBd)
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

%% define the optimization program
%   min_{d,\zeta_i}       Obj(d)
%          s.t.           Obj(d)=J(d)
%          s.t.           CVaR_w(alpha)<= 0
%          s.t.           d\in \Theta  
% where J(d) is a cost function  
% and w(d,\delta)=max_j g_j(d,\delta) is the worst-case performance fun.
%  CVaR_w(alpha) is estimated prescribing a Model for the uncertainty and
%  sampling realization of delta from it
Nd=length(dn); % get number of desing variables
LB=LBd; % lower bounds on d 
UB=UBd; %upper bounds on d  
X0=dn; % an initial guess 
obj= @(x) J_fun(x(1:Nd)); % The cost function  

%% model for the uncertainty
GMModel = fitgmdist(delta,size(delta,2)+3); % gaussian Mixture model 
 %% run fmincon optimzer 
% optionsGA = optimoptions('ga','PopulationSize',10000); % spq is needed to guaratnee high accuracy
%optionsGA
 % Scenario_Pf_constraint compute the constraints defined as w(d,delta_i)>\zeta_i 
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',50000,'Display','iter'); % spq is needed to guaratnee high accuracy
[X_opt,J_opt,exitflag,output]= fmincon(@(x) obj(x),... % objective function
             X0,[],[],[],[], LB, UB,...% initial guess and design bounds
              @(x) CVaR_Constraint(x,GMModel,w_fun,alpha,Nd),options); % non linear constraints and options
%[c]=CVaR_Constraint(X_opt,GMModel,w_fun,alpha,Nd); 
 %% Collect Results
Results.Jopt=J_fun(X_opt(1:Nd));
Results.dopt= X_opt(1:Nd); 
Results.Pf=mean(w_fun(X_opt,GMModel.random(1e6))>0);  
end