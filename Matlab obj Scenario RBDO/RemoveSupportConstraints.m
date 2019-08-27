 function   [Res]=RemoveSupportConstraints(obj,Support_Set,SRob,RemovalType)
            %  Support_Set= set of contraints to be removed
            %  e.g. Support_Set=SRob.Support_Set_Min;
            %  SRob= structure ouput of a ScenarioConstraints_ method 
            obj=obj.ScenarioSubset(); % restore original scenario set
            alpha=SRob.Opt.alpha;  Theta0=SRob.Opt.Theta0;  Tol=SRob.Opt.Tol; 
            if RemovalType==2 % if we remove all together
                disp('RemovalType  =2 (remove all together)')
                obj=obj.ScenarioSubset([setdiff(1:1:obj.Ndelta,Support_Set )]);
                [ScenarioRobustness_i]=obj.ScenarioConstraints_addMethod(Tol,alpha,Theta0);
                Res.ScenarioRobustness_i=ScenarioRobustness_i;
                Res.Theta_opt_i =ScenarioRobustness_i.Theta_opt;
                Res.epsilon_i =ScenarioRobustness_i.epsilon;
                Rel_Dopt_iremoved=obj.Compute_ReliabilityMetrics(Res.Theta_opt_i );
                Res.Pf_i =Rel_Dopt_iremoved.Pf;
                Res.Gmax_i=Rel_Dopt_iremoved.Gmax;
                Res.Gmax_K_i=max(Rel_Dopt_iremoved.G );
                Res.G_i=Rel_Dopt_iremoved.G;
            else  % if we remove one-at-a-time or sequentially
                if RemovalType==0 % one-at-a-time removal
                    disp('RemovalType  =0 (remove one-at-a-time)')
                elseif RemovalType==1
                    disp('RemovalType  =1 (remove sequentially)')
                end
                % prealocate for speed
                [ScenarioRobustness_i,Rel_Dopt_iremoved,G_i]=deal(cell(1,length(Support_Set)));
                [epsilon_iremoved,Pf_iremoved,Gmax_iremoved]=deal(zeros(1,length(Support_Set)));
                Gmax_K_iremoved=zeros(length(Support_Set),obj.Ng);
                Theta_newoptimal=zeros(obj.Nd+1,length(Support_Set));
                for i=1:length(Support_Set)
                    display(['Removed Suppport scenario: ' num2str(i) '/' num2str(length(Support_Set))])
                    if RemovalType==0 % one-at-a-time removal
                        obj=obj.ScenarioSubset(); % restore original scenario set
                        obj=obj.ScenarioSubset([setdiff(1:1:obj.Ndelta,Support_Set(i))]); % keep all but scenario i
                    elseif RemovalType==1 % Sequential removal
                        obj=obj.ScenarioSubset(); % restore original scenario set
                        obj=obj.ScenarioSubset([setdiff(1:1:obj.Ndelta,Support_Set(1:i))]); % keep all but scenario i
                    end
                    %[Theta_opt_lessi, Gamma_min_iremoved(i), exitflag(i), ~] = RBDO.Scenario_RBDO(alpha, Theta0, A, B, Aeq, Beq, LB, UB, options);
                    [ScenarioRobustness_i{i}]=obj.ScenarioConstraints_addMethod(Tol,alpha,Theta0);
                    Theta_opt_iremoved=ScenarioRobustness_i{i}.Theta_opt;
                    epsilon_iremoved(i)=ScenarioRobustness_i{i}.epsilon;
                    Rel_Dopt_iremoved{i}=obj.Compute_ReliabilityMetrics(Theta_opt_iremoved);
                    Theta_newoptimal(:,i)=Theta_opt_iremoved;
                    Pf_iremoved(i)=Rel_Dopt_iremoved{i}.Pf;
                    Gmax_iremoved(i)=Rel_Dopt_iremoved{i}.Gmax;
                    for k=1:obj.Ng
                        Gmax_K_iremoved(i,k)=max(Rel_Dopt_iremoved{i}.G(:,k));
                    end
                    G_i{i}=Rel_Dopt_iremoved{i}.G;
                end
                Res.ScenarioRobustness_i=ScenarioRobustness_i;
                Res.Theta_opt_i =Theta_newoptimal;
                Res.epsilon_i =epsilon_iremoved;
                Res.Pf_i =Pf_iremoved;
                Res.Gmax_i =Gmax_iremoved;
                Res.Gmax_K_i =Gmax_K_iremoved;
                Res.G_i=G_i;
            end
 end 