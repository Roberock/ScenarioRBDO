function [c,ceq]=Non_Linear_ScenarioProgramConstraints(obj,theta,alpha,Gexamined,Gdeltaidx)
%% This function evaluates the scenario constraints for programs SP1 and SP3
% alpha \in [0,1[ identifies the perentile of scenarios to keep
% theta=[d,gamma] d is the desing gamma is the upper bound on the worstcase
% obj= scenario RBDO object

Pcentile=(1-alpha)*100; % define percentile of samples to keep
ceq = [];  % no equality constraints ceq = 0
g = obj.Compute_Gfun(theta(1:obj.Nd));
gnorm=g./(1+abs(g));%% normalize g values in [-1 1]
Gamma=theta(obj.Nd+1:end);

if length(Gamma)==1 % Scenario program SP3
    W=max(gnorm,[],2);% the worst cases
    C_delta=(W<=prctile(W,Pcentile)); % select set of scenarios C_delta={delta: Fw (delta)<1-alpha}
    c_worstcase=max(W(C_delta));
    c=c_worstcase-Gamma;%% Now build inequality constraints
    %c=c*1e8; % harden soft non-linear constraints (check if feasibility in fmincon weights between obj minimization and constraints violation)
    
elseif length(Gamma)>1 % % Scenario program SP3
    c=zeros(obj.Ng,1);
    for ng=1:obj.Ng
        if any(Gexamined==ng) % this is needed for a selective remove of scenarios from g_j
            if isempty(Gdeltaidx) % if there are no idices worst case is safe g=0
                c_worstcase=-eps;
            else
                Gng=gnorm(Gdeltaidx,ng);% the worst cases for ng
                C_delta=(Gng<=prctile(Gng,Pcentile)); % select set of scenarios C_delta={delta: Fgng (delta)<1-alpha}
                c_worstcase=max(Gng(C_delta));
            end
            c(ng)=c_worstcase-Gamma(ng);%% Now build inequality constraints
        else
            Gng=gnorm(:,ng);% the worst cases for ng
            C_delta=(Gng<=prctile(Gng,Pcentile)); % select set of scenarios C_delta={delta: Fgng (delta)<1-alpha}
            c_worstcase=max(Gng(C_delta));
            c(ng)=c_worstcase-Gamma(ng);%% Now build inequality constraints
        end
    end
end


end