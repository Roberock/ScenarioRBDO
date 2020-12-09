function [CVaR]=CVaR(Samples,alpha)
%% Compute the CvaR by splitting probability atom
% CVar= \lambda_\alpha* VaR_\alpha+(1-\lambda_\alpha)E[X|X\geq VaR_\alpha]
% \lambda_\alpha= (F_X(VaR_alpha)-\alpha)/(1-\alpha)
% alpha \in [0,1] = confidence level
[F,X]=ecdf(Samples);
%XVaR=prctile(Samples,alpha );

i=find(F >= alpha  , 1, 'first');
XVaR=X(i);
if i==length(F) % if F_X(VaR_\alpha)=1   VaR_\alpha=max(Samples)
    CVaR=XVaR;
    % CVaR_minus=XVaR;
    % CVaR_plus=XVaR;
else % split the probability atom to assure continuity of CVaR for sample/discontinuous distribution
    lambda_alpha=( F(i)-alpha )/(1-alpha);
    %CVaR_minus=mean(X(X>=XVaR));
    CVaR_plus=mean(X(X>=XVaR));
    CVaR=lambda_alpha*XVaR+(1-lambda_alpha)*CVaR_plus;
end

end