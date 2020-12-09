function [c,ceq,gradc,gradceq]=CVaR_Constraint(x,GMModel,g_fun,alpha,Nd)
ceq=[]; rng 'default'
delta_GM=GMModel.random(1e5); % 10^5 samples from the GMixture model
G=g_fun(x(1:Nd),delta_GM); % reliability perform fun
W=max(G,[],2); % worst-case perform fun
[CVAR]=CVaR(W,alpha); % Conditional Value at Risk  
c=CVAR ; 
end