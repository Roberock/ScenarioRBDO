classdef ScenarioRBDO
    %% ScenarioRBDO: this class is a collection of methods and proprieties
    % and performs Reliability-based-desing-optimization using Scenario
    % Theory. The class is equipped with different optimization programs,
    % robustness assessment methods (to identify support sub-saples)
    % outlier removal methods and plotting capability.
    
    %%   Please refer to the following paper for further reading:
    %   Proceedings of the ASME 2019  Dynamic Systems and Control Conference, October 8-11, 2019, Park City, UT, USA
    %   RELIABILITY-BASED DESIGN BY SCENARIO OPTIMIZATION
    %   Cor. Author: Roberto Rocchetta, National Institute of Aerospace (NIA),  Hampton, VA, USA.
    %   Other Authors: Luis G. Crespo  andSean P. Kenny Dynamic Systems and Control Branch, NASA Langley Research Center, Hampton, VA, USA
    
    %%  %  %  %  %  % Example. How to use:%  %  %  %  % %  %  %  %  % %  %  %  %  %
    %          dn=[1,1,1]; %nominal desing
    %          delta=normrnd(0,1,[1e3,2]);  % 1e3 scenarios
    %          g_fun1=@(delta,d) +1/d(1)*(delta(:,2))+1/d(2)*(delta(:,1))-d(3);
    %          g_fun2=@(delta,d) +d(1)*(delta(:,1))-1/d(2)*(delta(:,2))-d(3);
    %          g_fun=@(delta,d) [g_fun1(delta,d) g_fun2(delta,d)]; % performance function
    
    %      %  %  %  % Now define the optimizer, bounds, constraints, and options %  %  %  %  %
    %         OptimizerData.LB=[.5 .5 .5  ,-inf]; %  constraints on [d,gamma]
    %         OptimizerData.UB=[2 2 2 ,+inf]; %  constraints on [d,gamma]
    %         OptimizerData.options= optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-6);
    %         [OptimizerData.A, OptimizerData.B, OptimizerData.Aeq, OptimizerData.Beq]=deal([]);
    %      %  %  %  %  %  Construct scenario reliability object %  %  %  %  %%  %  %  %  %
    %         RBDO=ScenarioRBDO('delta',delta,'dn',dn,'g_fun',g_fun,'OptimizerData',OptimizerData);
    
    %      %  %  % Analyse reliability of the nominal desing %  %  %  %  %  %  %  %
    %         Rel_dn=RBDO.Compute_ReliabilityMetrics(dn); display(Rel_dn)
    
    %      %  %  % Optimize desing (SP1 using all the scenarios) %  %  %  %  %  %  %
    %         alpha=0;  gamma0=0; Theta0=[dn gamma0];
    %         [Theta_opt_0, ~, ~, ~] = RBDO.Scenario_RBDO(alpha, Theta0);
    
    %      %  %  % Robustnes of the optimal desing %  %  %  %  %
    %         Tol=1e-4; % tollerance  the removed scenario
    %         [Robustness_opt0]=RBDO.ScenarioConstraints_addMethod(Tol,alpha,Theta0);
    
    %      %  %  % Optimize desing (SP1 removing 5 % of the scenarios)%  %  %  %  %
    %         alpha=0.05;  gamma0=0; Theta0=[dn gamma0];
    %         [Theta_opt_005, ~, ~, ~] = RBDO.Scenario_RBDO(alpha, Theta0)
    %  %   %  %  % %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
    %% START Proprieties
    properties
        d_opt     % [d^*] optimal desing leading to min gamma
        Sub_delta  %a subset of delta [Nsubsamples x N_scenario_dimension] Nsubsamples=<Nsamples
        OptimizerData % Structure with the following matlab optimizer input: A,B,Aeq,Beq,LB,UB,options
    end
    
    properties (SetAccess=protected)  % protected proprieties
        delta     % scenario matrix [Nsamples x N_scenario_dimension]
        g_fun     % performance function handle (g>0 defines failure) e.g g=@(delta,d) g_fun(delta,d) g is a [Nsamplex x Ng] matrix ;
        dn        % nomianl desing vector [1 x Ndesign_values]
    end
    
    properties (Dependent=true) % dependant proprieties
        Ng % number of performance functions
        Nd % number of design variables
        Ndelta % number of scenarios
        NSubdelta % number of scenarios
        Nrv % number of random variables (scenarios dimension)
    end
    
    methods
        
        %% CONSTRUCTOR
        function XSBRO   = ScenarioRBDO(varargin)  % method to define an object of type Scenario_Based_Reliability
            if nargin==0;  return % Create an empty object
            end
            for k=1:2:length(varargin)
                switch lower(varargin{k})
                    case {'delta','scenario','scenarios','samples'}
                        XSBRO.delta=varargin{k+1};
                        if size(XSBRO.delta,1)<size(XSBRO.delta,2)
                            warning('dimensionality > N of scenarios, you should input a [Nsamples x N_scenario_dimension] matrx')
                        end
                    case {'dn','desing','nomdes'};  XSBRO.dn=varargin{k+1};
                    case {'g_fun','performance','g'}; XSBRO.g_fun=varargin{k+1};
                    case {'sub_delta','subset','subsetscenario'}; XSBRO.Sub_delta=varargin{k+1};
                    case {'theta_opt','d_opt'}; XSBRO.d_opt=varargin{k+1};
                        if length(XSBRO.d_opt)~=length(XSBRO.dn); error('d_optimal differentl length than dnominal');end
                    case {'optimizerdata','optsetting'}; XSBRO.OptimizerData=varargin{k+1};
                    otherwise; error( 'PropertyName %s is not valid ', varargin{k});
                end
            end
            assert(~isempty(XSBRO.delta),'We need scenarios, delta can not be empty');
            assert(~isempty(XSBRO.g_fun),'We need a performance evaluation model, g_fun can not be empty');
            if isempty(XSBRO.Sub_delta)
                XSBRO=XSBRO.ScenarioSubset(); % if is empty, fill with the original scenario set
            end
        end     %of constructor
        %% EDIT OBJECT
        function obj=ScenarioSubset(obj,IDX,varargin)
            % IDX indices of the scenarios to keep in a subset
            if nargin == 2; obj.Sub_delta=obj.delta(IDX,:);
            elseif nargin == 1; obj.Sub_delta=obj.delta;
            end
        end
        %% GET DEPENDENT PROPRIETIES
        function Ndval = get.Nd(obj) % get number of performance functions
            if size(obj.dn,1)~=1;  error('Desing shoul be a vector [1 x Ndesign_variables)]');end
            Ndval=length(obj.dn); end
        function Ndeltaval = get.Ndelta(obj) % get number of performance functions
            Ndeltaval=size(obj.delta,1); end
        function Nrval= get.Nrv(obj) % get number of performance functions
            Nrval=size(obj.delta,2); end
        function Ngval = get.Ng(obj) % get number of performance functions
            Ngval=length(obj.g_fun(obj.delta(1,:),obj.dn)); end
        function NSubdeltaval = get.NSubdelta(obj) % get number of performance functions
            NSubdeltaval=size(obj.Sub_delta,1); end
        
        %% RELIABILITY ASSESSMENT METHODS
        function g = Compute_Gfun(obj,d)  % estimates Pf for a desing d and given the available scenario set delta
            g=obj.g_fun(obj.Sub_delta,d(1:obj.Nd));
        end
        function Pf = Compute_FailureProbability(obj,d) %Compute_FailureProbability estimates Pf for a desing d and given the available scenario set delta
            g=obj.g_fun(obj.Sub_delta,d(1:obj.Nd));
            Pf=mean(any(g>0,2));
        end
        function Gmax = Compute_maxG(obj,d)  %Compute_maxG computes the maximum g given the available scenario set delta
            g=obj.g_fun(obj.Sub_delta,d(1:obj.Nd));
            Gmax=max(max(g));
        end
        function W = Compute_W(obj,d) %Compute W(delta)= (max_J G_j(delta,d))
            g=obj.g_fun(obj.Sub_delta,d(1:obj.Nd));
            W=max(g,[],2);
        end
        function Reliability = Compute_ReliabilityMetrics(obj,d)
            %Reliability Summary of the reliability for desing d
            g=obj.g_fun(obj.Sub_delta,d(1:obj.Nd));
            Reliability.G=g; Reliability.Gmax=max(max(g));
            Reliability.Gmax_g=(max(g)); Reliability.Pf=mean(any(g>0,2));
            Reliability.W=max(g,[],2); Reliability.Pf_g=mean(g>0);
        end%Compute various the reliability metrics
        
        %%  OPTIMIZATION METHODS
        
        % NONLINEAR CONSTRAINTS METHOD
        function  [c,ceq]=Non_Lin_Cons_SP(obj,theta,alpha,Gexamined,Gdeltaidx,varargin)
            % this function evaluates scenario constraints for the programs SP1 and SP3
            if nargin<=3; Gexamined=[];Gdeltaidx=1:1:obj.NSubdelta;
            end
            [c,ceq]=Non_Linear_ScenarioProgramConstraints(obj,theta,alpha,Gexamined,Gdeltaidx);
        end
        
        % OPTIMIZERS
        % method (d*,gamma*)=SP_1(D_delta,alpha)
        function [Theta_opt,Gamma_min,exitflag,output] = Optimize_SP1(obj,alpha,Theta0)
            %  This method solve scenario program:
            %  (d*,\gamma*) =arg min{\gamma:  F_{w(d,\mathcal{D}_{\delta})}^{-1}(1-\alpha)  \leq \gamma \right \}
            
            if alpha<0 || alpha>=1; error('alpha should be 0=<alpha<1 ');end
            if length(Theta0)==(obj.Nd);  Theta0=[Theta0 +1];
            elseif length(Theta0)<(obj.Nd); warning('length guess desing < Number desing variables....guess set = dnominal'); Theta0=[obj.dn +1];  end
            
            % This function perfrom reliability based design optimization by using a scenario theoretical approach
            [A,B,Aeq,Beq, LBd, UBd, options]=getFminconSet(obj); % get optimization
            LB=[LBd,-1 ]; % add bounds on gamma
            UB=[UBd,+1 ]; % add bounds on gamma
            
            [Theta_opt,Gamma_min,exitflag,output]= fmincon(@(x) x(end),... % objective function
                Theta0,A,B,Aeq,Beq, LB, UB,...% initial guess and design bounds
                @(x) obj.Non_Lin_Cons_SP(x,alpha,[],[]), options); % non linear constraints and options
        end
        % method (d*,gamma*)=SP_2(D_delta)
        function [Theta_opt,Pf_min] = Optimize_SP2(obj)
            %  This method solve scenario program:
            %  (d*,\gamma*) =arg min{\gamma:  P_f \leq \gamma \right \}
            [~,~,~,~, LBd, UBd, ~]=getFminconSet(obj); % get optimization
            Pfminfun= @(x) mean(any(obj.Compute_Gfun(x)>0,2));
            gaDat.FieldD=[LBd;UBd];
            gaDat.Objfun=Pfminfun;
            % Execute GA
            gaDat=ga(gaDat); % gradient based methods are inapplicable, use a genetic algorithm
            % Result are in
            Theta_opt =gaDat.xmin;
            Pf_min =gaDat.fxmin; % this is the objective function Pf=gamma
        end
        % method (d*,gamma*)=SP_3(D_delta)
        function [Theta_opt,Gamma_min,exitflag,output] = Optimize_SP3(obj,alpha,Theta0,Gexamined,Gdeltaidx,varargin)
            %  This method solve scenario program:
            %  (d*,\gamma*) =arg {sum(max(\gamma_j,0)): F_{g_j(d,\mathcal{D}_{\delta})}^{-1}(1-\alpha)  \leq \gamma_j, j=1,..,ng \right \}
            if nargin<=4;Gdeltaidx=1:1:obj.NSubdelta; end % if not provided consider all the scenarios
            if nargin<=3;Gdeltaidx=1:1:obj.NSubdelta; Gexamined=[]; end
            if alpha<0 || alpha>=1; error('alpha should be 0=<alpha<1 ');end
            if length(Theta0)==(obj.Nd);  Theta0=[Theta0 ,ones(1,obj.Ng)];
            elseif length(Theta0)<(obj.Nd); warning('length guess desing < Number desing variables....guess set = dnominal'); Theta0=[obj.dn  ,ones(1,obj.Ng)];  end
            
            
            [A,B,Aeq,Beq, LBd, UBd, options]=getFminconSet(obj); % get optimization
            LB=[LBd,- ones(1,obj.Ng)]; % add bounds on gamma
            UB=[UBd,+ ones(1,obj.Ng)]; % add bounds on gamma
            
            Weights=ones(1,obj.Ng); % weights on the reliability performance functions g_j
            NormalizedWeights=Weights./sum(Weights); % normalize weights
            
            % run optimizer
            [Theta_opt,Gamma_min,exitflag,output]= fmincon(@(x) sum(max(NormalizedWeights.*x(obj.Nd+1:end),zeros(1,obj.Ng))),... % objective function
                Theta0,A,B,Aeq,Beq, LB, UB,...  % optimization inputs
                @(x) obj.Non_Lin_Cons_SP(x,alpha,Gexamined,Gdeltaidx), options);  % non-linear constraints and options
        end
        
        %% ROBUSTNESS ASSESSMENT METHODS
        % Method to get the robustness parameter epsilon given a confidence beta and number of supports k
        function  epsilon=getEpsilon(obj,k,beta)
            epsilon=getConfidence_nonconvex(k,beta,obj.Ndelta);
        end
        % ROBUSTNESS OF SP1
        % find support constraints using add Worst-Case-Scenario method
        function [Rbst]=ScenarioConstraints_addMethod(obj,Tol,alpha,Theta0)
            %this method identifies scenario constraints by adding the Worst
            % Case scenario i=arg_i max_j(g_j(Theta,delta(i)); for the candidate design Theta
            % It works as follows:
            % Phase 1 (add): The method converges to a reducible set of constraints
            % Phase 2 (remove): The support set is further reduced to a non-reducible support set
            % using the standard remove-one-at-a-time method
            disp('run Scenario RBDO');tic
            [Theta_opt,~,~,~] = obj.Optimize_SP1(alpha,Theta0); %First we Optimize the desing with SP1
            CompOptimization=toc;% computational time for the first optimization
            disp('Identify scenario constraints by adding Worst-Case-Scenario one-at-a-time')
            obj.OptimizerData.options.Display='off'; %obj.OptimizerData.options.MaxIterations=200;
            Rbst= Support_Constraints_addMethod_SP1(obj,Theta_opt,alpha,Theta0,Tol); % run phase 1 and phase 2
            % add info to the results  strucure
            Rbst.Theta_opt=Theta_opt;
            Rbst.Opt.Theta0=Theta0; Rbst.Opt.Tol=Tol; Rbst.Opt.alpha=alpha;
            Rbst.beta=1e-8;% default confidence
            Rbst.epsilon=getEpsilon(obj,Rbst.Cardinality,Rbst.beta);%get reliability index
            Rbst.ComputationalTime.Optimization=CompOptimization;
        end
        % find support constraints using remove one-at-a-time
        function [Rbst]=ScenarioConstraints_removeMethod(obj,Tol,alpha,Theta0)
            %this method removes scenario one by one, it checks weather scenario support constraints
            % by comparing the optimal reference desing Theta_opt with the solution obtained removing the scenario
            disp('running Scenario RBDO');tic
            [Theta_opt,~,~,~] = obj.Optimize_SP1(alpha,Theta0);% First optimize
            CompOptimization=toc;% computational time for the first optimization
            disp('Identify scenario constraints using one-at-a-time scenario removal methods')
            obj.OptimizerData.options.Display='off';  % hide optimization display
            [Rbst]=Support_Constraints_removeMethod_SP1(obj,alpha,Theta_opt,Theta0,Tol); % run the method
            % add info to the results  strucure
            Rbst.Theta_opt=Theta_opt;
            Rbst.Opt.Theta0=Theta0; Rbst.Opt.Tol=Tol; Rbst.Opt.alpha=alpha;
            Rbst.beta=1e-8; % default confidence
            Rbst.epsilon=getEpsilon(obj,Rbst.Cardinality,Rbst.beta); %get reliability index
            Rbst.CompTimeOptimizer=CompOptimization;
        end
        % ROBUSTNESS OF SP3
        function [Rbst]=ScenarioConstraintsSp3_removeMethod(obj,G2examin,Tol,alpha,Theta0,Gdeltaidx)
            % Tol= tollerance e.g. 1e-6
            % Theta0= inital guess
            % G2examin= the performance function index to examin, e.g. G2examin=2, we examin g2 in g=[g1 g2 g3]
            % alpha=0 keeps all the scenarios ( alpha>0 momentarely no theory to assess the program robustness)
            
            %  First RUN OPTIMIZER
            % Scenario_RBDO  This method run a Reliability-based Desing-Optimization with scenario theoryalphaalpha
            % Solve scenario program: {\gamma:  F_{w(d,\mathcal{D}_{\delta})}^{-1}(1-\alpha)  \leq \gamma \right \}
            [Theta_opt,~,~,~]=obj.Optimize_SP3(alpha,Theta0,G2examin,Gdeltaidx);
            [Rbst]=Support_Constraints_removeMethod_SP3(obj,Theta_opt,alpha,Theta0,Tol,Gdeltaidx,G2examin);Rbst.Opt.Theta0=Theta0; Rbst.Opt.Tol=Tol; Rbst.Opt.alpha=alpha;
            Rbst.Theta_opt=Theta_opt;
            Rbst.G2examin=G2examin;
            Rbst.beta=1e-8;% default confidence
            Rbst.epsilon=getEpsilon(obj,Rbst.Cardinality,Rbst.beta);%get reliability index
        end 
        function [Rbst]=ScenarioConstraintsSp3_addMethod(obj,G2examin,Tol,alpha,Theta0)
            Gdeltaidx= (1:1:obj.NSubdelta)';
            [Theta_opt,~,~,~]=obj.Optimize_SP3(alpha,Theta0,G2examin,Gdeltaidx);
            [Rbst]=Support_Constraints_addMethod_SP3(obj,Theta_opt,alpha,Theta0,Tol,G2examin);
            Rbst.Theta_opt=Theta_opt; 
            Rbst.beta=1e-8;% default confidence
            Rbst.G2examin=G2examin;
            Rbst.epsilon=getEpsilon(obj,Rbst.Cardinality,Rbst.beta);%get reliability index
        end
        
        %% OUTLIERS REMOVAL METHODS
        function [Res]=RemoveConstraints(obj,Support_Set,SRob,RemovalType,varargin)
            obj.OptimizerData.options.Display='off';
            if nargin<=3;RemovalType=2;
                warning('RemovalType set to default 2. input RemovalType=0 (remove one-at-a-time) =1 (remove sequentially) or =2 (remove all together)')
            end
            [Res]=RemoveSupportConstraints(obj,Support_Set,SRob,RemovalType);
        end
         
        %% PLOT METHODS
        function plot_deta_vs_G(obj,d,Idx_delta,Idx_G)
            % plot_deta_vs_G(Idx_delta,Idx_G)
            %Plot_ScenarioConstraints: the desing to plot
            % Idx_delta: 2-D the uncertainty dimensions to plot, e.g. Idx_delta=[1,2] 
            % Idx_G: 2  list of performance dimensions to plot e.g. Idx_G=[1] or more Idx_G=[1,2,5]
            if length(Idx_delta)~=2
                error('plot_deta_vs_G(d,Idx_delta,Idx_G) Idx_delta should be a vector with 2 scenario indices e.g. Idx_delta=[1,2]');  end
            if max(Idx_G)>obj.Ng || isempty(Idx_G)
                error('plot_deta_vs_G(Idx_delta,Idx_G) Idx_G is empty or exceeds the max number of g performances');  end
            
            figure
            g=obj.g_fun(obj.Sub_delta,d(1:obj.Nd));
            for i=1:length(Idx_G)
                subplot(1,length(Idx_G),i)
                scatter3(obj.Sub_delta(g(:,Idx_G(i))>0,Idx_delta(1)),obj.Sub_delta(g(:,Idx_G(i))>0,Idx_delta(2)),g(g(:,Idx_G(i))>0,Idx_G(i)),'r') %plot success cases
                hold on
                scatter3(obj.Sub_delta(g(:,Idx_G(i))<=0,Idx_delta(1)),obj.Sub_delta(g(:,Idx_G(i))<=0,Idx_delta(2)),g(g(:,Idx_G(i))<=0,Idx_G(i)),'*b')  %plot failure cases
                xlabel(['\delta_' num2str(Idx_delta(1))]);ylabel(['\delta_' num2str(Idx_delta(2))]);zlabel([['g_' num2str(Idx_G(i))] '(\delta,d)'])
                set(gca,'FontSize',20); grid on; box on
            end
        end
        
        function plot_detaIndex_vs_G(obj,d,IDX)
            % plot the the performance function vs the scenario realization
            % sorted by g_{IDX}
            
            % d: the desing vector
            % IDX: the index of the performance function (g_{IDX}) that we
            % would like to sort
             
            g=obj.g_fun(obj.Sub_delta,d(1:obj.Nd));
            [~,deltaIdx]=sort(g(:,IDX));
            G=g(deltaIdx,:); %sort (generally provides a better graphical output)
            LEGEND={};STYLES={'-b','--r','.-k','-r','--k','.-b','-k','--b','.-r'};
            for k=1:obj.Ng
                hold on
                plot(G(:,k),STYLES{randi(length(STYLES))})
                LEGEND{k}= ['g_' num2str(k)];
            end
            grid on
            xlabel(['Scenario ID sorted by g_' num2str(IDX)]);ylabel('g_i');
            set(gca,'FontSize',20); grid on; box on
            legend(LEGEND)
        end
        
        function plot_2D_SafeFailDomains_and_Scenarios(obj,d,DeltaIndex,GIndex,varargin)% Output plot: a 2d vew of the different desing and performance functions
            % d=desing; DeltaIndex= 2 scenario indices; GIndex= g indixes
            % plot failure-safe region
            % Input:  2-D dimension of delta, 1 or more designs, 1 or more performance function idices
            
            if max(DeltaIndex)>obj.Nrv
                error(['DeltaIndex is a 2-d vector with scenario indices, max idx should be less than ' num2str(obj.Nrv)])
            end
            if nargin==4; PlotGall=0;
            elseif nargin==3; PlotGall=1;GIndex=[] ;
            elseif nargin==2; PlotGall=0; DeltaIndex=[1,2];GIndex=1:1:obj.Ng;
            elseif nargin==1; PlotGall=0; d=obj.dn; DeltaIndex=[1,2];GIndex=1:1:obj.Ng;
            end
            
            plot_2D_SafeFailDomain(obj,d,DeltaIndex,GIndex,PlotGall)
        end
        
        function Plot_ScenarioConstraints(obj,d,deltaIdx,varargin)% plot scenario constraints in the max(g)-d space
            if nargin<3;  disp(' scenario indices not provided: deltaIdx includes all the scenarios')
                deltaIdx=1:1:obj.NSubdelta;
            elseif nargin<2;  disp('desing and scenario indices not provided: d was set to the nominal value dn and deltaIdx includes all the scenarios')
                d=obj.dn;  deltaIdx=1:1:obj.NSubdelta;
            end
            Plot_ScenarioConstraints_dvsw(obj,d,deltaIdx)
        end
    end
end

%%%%%% THE LIST OF METHODS

%  Reliability Methods: INPUT: a system desing (d)

%      Compute_Gfun(d): evaluates the performance function g=[g1,..,gNg]
%      Compute_FailureProbability(d): evaluates the overall Pf on the available scenarios
%      Compute_maxG(d):
%      Compute_W(d)
%      Compute_ReliabilityMetrics(design)

%  Scenario Optimization Methods:
 
%           Th program SP2 minimize Pf (to this end it does not induce constraints)
%   NLCon-  NonLinearConstraint(theta,alpha,Gexamined,Gdeltaidx) % evaluate
%           non linear constraints for program SP1 or SP3

%  SP-1)    Optimize_SP1(alpha,Theta0):  minimizes alpha percentile of w (using fmincon)
%  SP-2)    Optimize_SP2(): minimizes Pf given-data (using GA);
%  SP-3)    Optimize_SP3(alpha,Theta0,Gexamined,Gdeltaidx):minimizes the sum
%           of the alpha percentiles of each requitement g_j with j=1,..,Ng  (using fmincon)

% Robustness method
%   getEpsilon(k,beta)  gets the non-convex robustness given cardinality k  and confidence beta (scenario size N is within the object)
 
%   ScenarioConstraints_addMethod (for SP1 and SP3)
%   ScenarioConstraints_removeMethod (for SP1 and SP3)

%  Outlier Removal Method
%   RemoveConstraints re-optimizes removing a list of scenario from the initial data set

%  Plots
%    plot_deta_vs_G scatter: the scenarios in the uncertainty space (2-D) vs the performance function realizations
%    plot_detaIndex_vs_G sort and plot the gj and the scenarios indices w.r.t. one of the reliability requirement
%    plot_2D_SafeFailDomains_and_Scenarios % plot failure and safe regions
%    Plot_ScenarioConstraints : plot a list of scenario constraints 