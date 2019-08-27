function  [Rbst]=Support_Constraints_removeMethod_SP3(obj,Theta_opt,alpha,Theta0,Tol,Gdeltaidx,G2ex)

CountIter=0; 
Support_Set=[];
RMSE=zeros(1,obj.NSubdelta);

while true
    ToDo=setdiff(Gdeltaidx,Support_Set);
    if isempty(ToDo); break % did we check all the scenarios for gj? yes. Then, break while
    else % No. Then, extract randomly a new constraint sceanrio to be tested
        idx_i=randi(length(ToDo));
%         while any(idx_i==Support_Set) % it avoids repeating analysis on known support scenarios
%             idx_i=randi(length(Gdeltaidx));
%         end
        Scenario_i=ToDo(idx_i); % pick on of the remaining randomly
        Gdeltaidx(Gdeltaidx==Scenario_i)=[];
    end
    %% 2) run the optimization with reduced set of scenarios Gdeltaidx and on the requirement G2examin
    % Scenario_RBDO  This method run a Reliability-based Desing-Optimization with scenario theoryalphaalpha
    % Solve scenario program: {\gamma:  F_{w(d,\mathcal{D}_{\delta})}^{-1}(1-\alpha)  \leq \gamma \right \}
    [Theta_removed_i,~,~,~]=obj.Optimize_SP3(alpha,Theta0,G2ex,Gdeltaidx);
    %% 3) Check Errors w.r.t. the reference optimized solution Theta_opt
    Error=sqrt((Theta_opt-Theta_removed_i).^2)./Theta_opt;
    CountIter=CountIter+1; % conunt iterations
    % RMSE(CountIter)= sqrt(mean((Theta_removed_i- Theta_opt).^2));
    RMSE(CountIter)= sqrt(mean((Theta_removed_i([1:obj.Nd,G2ex])- Theta_opt([1:obj.Nd,G2ex])).^2));
    %ISequalDesign=all(Error<Tol);% logic 1, then the error are small enough
    ISequalDesign=RMSE(CountIter)<Tol;
    if ISequalDesign~=1 % the removed scenario changes the solution, it is a support constraint
        Gdeltaidx=[Gdeltaidx  Scenario_i]; % put the scenario back
        Support_Set=union(Support_Set,Scenario_i); % include the scenario in the list of support constraints
    end
end
Rbst.Support_Set=Support_Set;
Rbst.Cardinality=length(Support_Set);
Rbst.CountIter=CountIter;
Rbst.ComputationalTime=toc;
Rbst.RMSE=RMSE(1:CountIter);

end