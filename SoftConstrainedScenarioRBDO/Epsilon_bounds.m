function epsilon=Epsilon_bounds(SN,N,beta,OptType,Nd)
%% Function to compute the violation probability bounds for different types of Scenario optimization programs
%
% Input:
% SN = number of support scenarios
% N = number of scenarios
% OptType = type of scenario bound
% OptType \in {'convex_priori','convex_w','convex_soft','convex_discard','nnconvex'}
% Nd = number of design variables
% Output:
% epsilon= the bound on the violation probabilty

if strcmp(OptType,'convex_priori')    % Convex programs: Nd<Nscenarios (a-priori beta dominated)
    
    
elseif strcmp(OptType,'convex_w&j')    % Convex programs: (a-posteriori Wait & Judge)
    
    out = getWaitandJudgeEpsilon(SN,N,beta);
    epsilon=out(end);
elseif strcmp(OptType,'convex_soft')    % Convex programs: (a-posteriori with soft constraints)
    
    % Reference Article
    % Title: 'Risk and complexity in scenario optimization'
    % Authors S. Garatti and ·M. C. Campi  Journal Mathematical Programming
    % DOI https://doi.org/10.1007/s10107-019-01446-4
    [epsL, epsU] = epsLU(SN,N,beta);
    epsilon=[epsL,epsU];
    
elseif strcmp(OptType,'convex_discard')    % Convex programs: (a-posteriori with removed constraints SN)
    EPS=linspace(0.01,0.99,1000);
    beta=getConfidence_ConvexDiscard(N,EPS,SN,Nd);
    epsilon=EPS(find(beta<=10^-8,1,'first'));
    
elseif strcmp(OptType,'nnconvex')
    % Non-Convex program: (a-posteriori)
    epsilon=getConfidence_nonconvex(SN,N,beta);
end
end