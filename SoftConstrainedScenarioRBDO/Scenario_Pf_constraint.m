function [c,ceq]=Scenario_Pf_constraint(x,delta,performance_fun,Nd)

%% Constraint function for Scenario program with slack variables (softened constraints)
% x= [design vector ,  vector of slack variables \zeta]  [1 x Nd+N]
% delta= matrix with the samples from the uncertainty [N x Ndelta] 
% performance_fun= a function f(x,delta)
% Nd= number of desing variables 

% thes are the orignial constraints f(x,delta^{(i)}) <= 0 i=1,...,N
% thes are the softened constraints f(x,delta^{(i)}) <= \zeta^{(i)}, i=1,...,N
ceq=[];
% evaluate f(x,delta)
G=performance_fun(x(1:Nd),delta);
SlakVariables=reshape(x(1+Nd:end),size(G,1),size(G,2));
% apply slack variables softening
SlackenG=G-SlakVariables;
c=SlackenG(:);
 
end