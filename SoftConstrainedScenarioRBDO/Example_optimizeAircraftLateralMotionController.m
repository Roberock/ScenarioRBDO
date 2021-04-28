clc; 
clear variables;
close all
%% Example of an aircraft lateral motion controller  
addpath([pwd '\G_functions']) % add path with the reliability functions (system-solvers)

%% Example
% dn=[0.3075 -0.3164 -0.0973 -0.0188 -0.3164 0.5822 -0.0703 -0.0993 ...
%     -0.0973 -0.0703 0.2277 0.2661 -0.0188 -0.0993 0.2661 0.7100 -0.0191 ...
%     0.2733 -0.0920 0.4325 0.0803 -0.3821 0.4496 -0.2032]; % nominal design from previous studies
dnominal=[0.3075 -0.3164 -0.0973 -0.0188 0.5822 -0.0703 -0.0993  0.2277 0.2661  0.7100... % P matrix components
    -0.0191 0.2733 -0.0920 0.4325 0.0803 -0.3821 0.4496 -0.2032];% W matrix components
Nd=length(dnominal); 
Ng=1; % number of reliability requirements
LBd=dnominal-abs(dnominal)*0.3;% lower bounds on the design parameters
UBd=dnominal+abs(dnominal)*0.3;% upper bounds on the design parameters

%% Define reliability, cost and worst-case reliability functions
g_fun = @(d,delta)    G_fun_lateralmotioncontroller(d,delta);  % it is the max eigneval of quadratic performance function
J_fun = @(d)          trace(reshape(d(1:16),4,4));  % Trace of the matrix P defines the cost
w_fun = @(d,delta)    max(g_fun(d,delta),[],2); % Worst-Case-Performance function

%% Define the data generating mechanism (unknown to the designer)
% for illustration we use normal distirbuted nose
% the family of the probabiiy distirbution does not matter....only the available samples are used within the optimizaiton.
N=150; %number of samples available for the optimization
Ntest=10^4;  % number of samples for testing the validity of the results

% nominal value of 13 uncertain factors
delta_mu=[ -2.93 -4.75 0.78 0.086 -0.11 0.1 -0.042 2.601 -0.29 -3.91 0.035 -2.5335 0.31]; 
% define a function for the Data-Generating-Mechanism 
DGM= @(N) normrnd(0,0.2 ,[N,13]).*delta_mu +delta_mu; 
 
%  generate samples 
delta=DGM(N); % scenarios for the RBDO 
delta_testing=DGM(Ntest); % scenarios for testing the results

%% Let us now compute the `true' (testing) reliability of the nominal desing
Rel_nominal_controller=ComputeReliabilityPerformance(dnominal, delta_testing,g_fun);  % compute the reliability for the testing data (as a reference)
disp(Rel_nominal_controller) % print statistics

%% RUN the SCENARIO RBDO 
VaR=0; % slect a value-at-risk level (lambda in the paper)
RHO=10^5; % cost associated to constraint violations
tic % computational time
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',5000,...
                    'display','off','UseParallel',false);
[Results_controller,~,~,~,optresults]=...
            ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dnominal,LBd,UBd,options);
CompTime=toc;  % save computational time 
% compute the reliability of the otpimized design using the testing data
Rel_sp_opt=ComputeReliabilityPerformance(Results_controller.dopt ,delta_testing,g_fun); 
sN=Results_controller.Support.Size; % get size of the support sceanarios
disp(Rel_sp_opt) % print statistics for the optimized design

% Robustness
beta=10^-4; % the confidence parameter (confidence level is given by 1-beta)
OutWeJ = getWaitandJudgeEpsilon(sN ,N,beta);% get Wait-and-judge reliability
epsilon_wait_and_judge =OutWeJ(end); % get Wait-and-judge reliability for the computed support size
[EpsilonLU(1), EpsilonLU(2)]= epsLU(sN,N, beta); % compute lower and upper bounds on the probability of exceeding the selected value at risk
 

disp(['True Failure Probability for the optimized design: Pf = '...
        num2str(Rel_sp_opt.Pf_all)])
disp(['Number of support constraints in the scenario RBDO problem '...
        num2str(Results_controller.Support.Size)]) 
 disp('Scenario bounds on the probabilty of exceeding the selected value-at-risk level') 
 disp(['For a VaR=0 the scenario-based bounds on the failure probability are ['...
        num2str(EpsilonLU) ']'])  
disp(['Obtained For a confidence level 1-beta = '...
        num2str(1-beta) ' and number of scenario constraints N = ' num2str(N) ])     

    
    %% RUN the SCENARIO RBDO 
VaR=-.5; % do we want to be over conservative? VaR=-0.5  (extra .5 on the safety margin)
RHO=10^5; % cost associated to constraint violations
tic % computational time
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',5000,'display','off','UseParallel',false);
[Results_controller_var2,X_opt,~,~,~]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dnominal,LBd,UBd,options);
CompTime=toc;  % save computational time
% compute the reliability of the otpimized design using the testing data
Rel_sp_opt_var2=ComputeReliabilityPerformance(Results_controller_var2.dopt ,delta_testing,g_fun); 
sN=Results_controller_var2.Support.Size; % get size of the support sceanarios


figure(1)
ecdf(Rel_nominal_controller.w); hold on;
ecdf(Rel_sp_opt.w)
ecdf(Rel_sp_opt_var2.w)

    %% RUN the SCENARIO RBDO 
VaR=-.3; % now VaR=-0.1 (extra .51on the safety margin)
RHO=10^5; % cost associated to constraint violations
tic % computational time
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',5000,'display','off','UseParallel',false);
[Results_controller_var3,X_opt,~,~,~]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dnominal,LBd,UBd,options);
CompTime=toc;  % save computational time
% compute the reliability of the otpimized design using the testing data
Rel_sp_opt_var2=ComputeReliabilityPerformance(Results_controller_var3.dopt ,delta_testing,g_fun); 
sN=Results_controller_var3.Support.Size; % get size of the support sceanarios



figure(2) 
ecdf(Results_controller_var2.Zopt);  hold on;
ecdf(Results_controller_var3.Zopt); 