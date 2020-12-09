function Epsilon=Get_Epsilon_Individual_Requirements(g_fun,delta,d_opt,beta,VaR)
%% This script computes the lower and upper bounds on the probability of violation for individual requirements
% Input
% beta = confidence parameter e.g. beta=10^-6
% d_opt = optimal desing
% delta = matrix of random samples
% g_fun = reliability performance function g_fun= @(d,delta) .....
% VaR= value-at-risk (e.g VaR=0 if we look at the probabilisty of failure)

% output
% Epsilon = lower and upper bounds on the violation probabilty for each requirement g_j with j=1,...,Ng

% Prob[ epsilon_lower <= V(d_opt) <= epsilon_upper} >= 1-beta
% V(d) =Prob[delta: g_j(d_opt,delta)>0]
%%
g=g_fun(d_opt,delta);
N=size(delta,1);
Ng=size(g,2);
SN_ind=sum(g>=VaR); % number of support scenarios for each g_j
Epsilon=zeros(Ng,2);
for i=1:Ng 
    [EpsilonLU(1), EpsilonLU(2)]= epsLU(SN_ind(i),N, beta);
    Epsilon(i,:)=[EpsilonLU(1), EpsilonLU(2)]';
end

end