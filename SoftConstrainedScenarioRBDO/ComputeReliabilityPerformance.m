function Rel=ComputeReliabilityPerformance(d,delta,g_fun)
Rel.g= g_fun(d,delta); % Worst-Case-Performances Across Samples (normalized [-1 1])
Rel.w= max(g_fun(d,delta),[],2); % Worst-Case-Performances Across Samples (normalized [-1 1])
%%% Useful Functions %% %% %% % %%
g_norm=@(d,delta) g_fun(d,delta)./(1+abs(g_fun(d,delta)));
Rel.g_normalized= g_norm(d,delta); % Performance Functions (normalized [-1 1])
Rel.w_normalized= max(g_norm(d,delta),[],2); % Worst-Case-Performances Across Samples (normalized [-1 1])
Rel.Pf_individual= mean(g_fun(d,delta)>=0); %  lambda-relaxed Failure Probability for individual requiremets
Rel.Pf_all= mean(max(g_fun(d,delta)>=0,[],2)) ; % lambda-relaxed Overall Failure Probability
Rel.g_max= max(Rel.g); % worst case for individual g for the desing d
Rel.w_max= max(Rel.w); % worst case scenario w for the desing d
Rel.w_p95= prctile(Rel.w,95); %  value at risk of w at alpha=0.95
Rel.CVAR_w95= CVaR(Rel.w,0.95); % conditional value at risk of w  at alpha=0.95
Rel.CVAR_alphaPf= CVaR(Rel.w,1-Rel.Pf_all); % conditional value at risk of w  at alpha=0.95
end