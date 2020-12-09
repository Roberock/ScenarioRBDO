function [c,ceq,gradc,gradceq]=CVaR_Constraint_Rockafellar2010(x,GMModel,g_fun,alpha,Nd,Nsamples)
% R.T. Rockafellar et al, "On buffered failure probability in design and optimization of structures
% Reliability Engineering and System Safety 95 (2010) 499–510


% min_{d,\zeta_0,\zeta} J(d) 
% s.t. \zeta_0+1/(N*alpha0)sum_{i=1}^N  \zeta^{(j)} \leq 0
% w(d,\delta^{i})-\zeta_0 \leq \zeta^{(j)}
% \zeta^{(j)} \geq 0
d=x(1:Nd);
zeta0=x(Nd+1);
zetaj=x(Nd+2:end);

ceq=[]; rng 'default'
%Nsamples=1e5;
delta_GM=GMModel.random(Nsamples); % 10^5 samples from the GMixture model
G=g_fun(d,delta_GM); % reliability perform fun
W=max(G,[],2); % worst-case perform fun
%alpha0=mean(W<0); % 1-Pf
alpha0=alpha;
c1= zeta0 + 1/(Nsamples*(1-alpha0))*sum(zetaj); 
c2= W-zeta0-zetaj';
c=[c1;c2];
end