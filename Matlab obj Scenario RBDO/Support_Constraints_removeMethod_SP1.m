function [ScenarioRobustness]=Support_Constraints_removeMethod_SP1(obj,alpha,Theta_opt,Theta0,Tol)
% initialize, preallocate
CountIter=0; % itertions counter
Scenarios=obj.Sub_delta; % scenario vectors
ScenarioRemaining=(1:1:size(Scenarios,1))'; % scenario idx to check
Checked=[]; % scenario idx checked
L_removed_i= [Scenarios ScenarioRemaining];
Support_Set=[]; % idx with support scenarios
RMSE=zeros(1,obj.Ndelta);
tic
while true
    ToDo=setdiff(ScenarioRemaining,Support_Set);
    if isempty(ToDo)% did we check all the scenarios? yes. Then, break while
        break
    else % No. Then, extract a new scenario to be checked
        idx_i=randi(length(ToDo));
        Scenario_i=ToDo(idx_i); % pick on of the remaining randomly
        Checked=union(Checked,Scenario_i);
    end
    L_removed_i(idx_i,:)=[]; % remove the scenario from the list 
    
    %% 2) run the optimization with reduced set of scenarios
    obj.Sub_delta=L_removed_i(:,1:end-1);
    Theta_removed_i= obj.Optimize_SP1(alpha, Theta0);
    %% 3) Check Errors
    Error=sqrt((Theta_opt-Theta_removed_i).^2)./Theta_opt;
    CountIter=CountIter+1;
    RMSE(CountIter)= sqrt(mean((Theta_removed_i- Theta_opt).^2));
    
    ISequalDesign=all(Error<Tol);
    if ISequalDesign~=1 % the removed scenario changes the solution, it is a support constraint
        L_removed_i=[L_removed_i; Scenarios(Scenario_i,:) Scenario_i]; % put the scenario back
        Support_Set=union(Support_Set,Scenario_i); % include the scenario in the list of support constraints
    end
    ScenarioRemaining=L_removed_i(:,end); % last column contains the scenario indices 
    %ToDo=setdiff(Remaining_Scenarios,Examined_Scenario_Set);
      
end
ScenarioRobustness.Support_Set=Support_Set;
ScenarioRobustness.Cardinality=length(Support_Set);
ScenarioRobustness.CountIter=CountIter;
ScenarioRobustness.ComputationalTime=toc;
ScenarioRobustness.RMSE=RMSE(1:CountIter);
end