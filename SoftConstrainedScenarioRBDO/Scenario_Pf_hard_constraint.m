function [c,ceq]=Scenario_Pf_hard_constraint(x,delta,performance_fun)

%% Constraint function for Scenario program with slack variables (softened constraints)
% x= [design vector ,  vector of slack variables \zeta]  [1 x Nd+N]
% delta= matrix with the samples from the uncertainty [N x Ndelta] 
% performance_fun= a function f(x,delta)
% Nd= number of desing variables 
 
ceq=[]; 
G=performance_fun(x,delta);  
c=G(:);
 
end