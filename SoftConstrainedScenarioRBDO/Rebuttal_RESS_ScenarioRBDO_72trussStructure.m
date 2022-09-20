%% REBUTTAL: Scenario-based RBDO optimization of the 72 truss, 4 floors story building

clc; clear variables; close all
%% Example of truss strure RBDO
addpath([pwd '\G_functions'])
addpath([pwd '\G_functions\Truss'])
N=500; % number of available samples
Ntest=10^4;% more samples for testing (these are not available in practice)
%% Example of 72-bar four level skeletal tower
% 1 [lb] = 0.45359237 [kg]
% 1 [in] = 2.54 [cm] =0.0254 [m]
% 1 [in^2] = 6.4516 [cm^2] =645.16 [mm^2]
% 1 [mm^2] = 0.00155 [in^2]
D=Data72; % Define the truss data (geometry, connectivity, etc)
LBd=D.LB;
UBd=D.UB+1;
%% START
Nd=16; % 16-desgin variables
%Ng=12;
RandomModel=2; % select model for the uncertianty
if RandomModel==1 % uncertainty in the 72 E modulus and in 12 (x,y,z) components of the loads in nodes 1-4
    DGM= @(N)  [unifrnd(-5*10^3,5*10^3,[N,12]) normrnd(0,10^5,[N,72])];
elseif RandomModel==2 % uncertainty in the 72 E modulus and in 6 components of the loads in nodes 1-4
    DGM= @(N)  [normrnd(5*10^3,500,[N,6]) normrnd(0,10^5,[N,72])];
elseif RandomModel==3 % uncertainty in the load at node 1 and in the Area...
    DGM= @(N) [normrnd(5*10^3,50,[N,3]) unifrnd(-0.01,0.01,[N,72])];
elseif RandomModel==4 % uncertainty in six loads, E modulus and Area...
    DGM= @(N) [normrnd(5*10^3,10,[N,6]) normrnd(0,10^3,[N,72]) unifrnd(-0.01,0.01,[N,72])];
end
g_fun = @(d,delta) g_fun_ST72(d,delta,D,RandomModel);
J_fun = @(d) J_fun_ST72(d,D);
w_fun = @(d,delta) max(g_fun(d,delta),[],2); % Worst-Case-Performance function

% data
delta=DGM(N);
delta_testing=DGM(Ntest);
%% a list of deterministic baseline designs (for the 72-truss)
% dn=[1.88, 0.52, 0.10, 0.10, 1.28, 0.52, 0.10, 0.10,  0.54, 0.52, 0.10, 0.10, 0.16, 0.55, 0.41, 0.58]; % baseline deterministic design
% dn1=[0.1, 0.682, 0.533, 0.661, 0.727, 0.67, 0.10, 0.10,  1.665, 0.668, 0.10, 0.10, 2.447, 0.667, 0.1, 0.1]; % https://doi.org/10.1016/j.advengsoft.2015.11.001  SORA-ICDE
% dn2=[0.1, 0.9, 0.5 , 0.7 , 0.9 , 0.9 , 0.3 , 0.3 ,  1.9 , 0.7, 0.3 , 0.10, 2.1, 0.5, 0.1, 0.1]; % see DOI 10.1007/s00366-015-0427-9
% dn3=[2.025, 0.533, 0.1, 0.1, 1.1567,0.5689, 0.10, 0.10,  0.5137, 0.479, 0.10, 0.10, 0.1579, 0.5501, 0.3449, 0.4984]; %  GP Adeli andKamal https://www.sciencedirect.com/science/article/pii/0045794986903275 \cite{ADELI1986501}
% dn_PSObest=[0.1615, 0.509, 0.496, 0.561, 0.514, 0.546, 0.1, 0.1095, 1.307, 0.519, 0.1, 0.1 , 1.742, 0.5185, 0.1, 0.1]; % PSO best Perez andBehdinan https://doi.org/10.1016/j.compstruc.2006.10.013
% dn_PSO=[1.7427, 0.518, 0.1,  0.1, 1.3079, 0.5193, 0.1, 0.1, 0.5142, 0.5464, 0.1, 0.1095 , 0.1615, 0.502, 0.4967, 0.5619]; % PSO https://reader.elsevier.com/reader/sd/pii/S004579491100277X?token=199A562FF09414EEAD6E9B34982AF2D86F1A3BC9669F5EC03B22CB71D2AEACC97D4610D167970FC19757F60ED6FCAD91
% % Designs reviewed in TABLE 7 of    https://doi.org/10.1016/j.advengsoft.2015.11.001
% %Hybrid Particle Swarm Optimization
% d_HSPO=[1.857, 0.505, 0.1, 0.1, 1.255, 0.503, 0.1, 0.1, ...
%     0.496, 0.506, 0.1, 0.1, 0.1, 0.524, 0.4, 0.534];
% %Harmony Search;
% d_HS=[1.79, 0.521, 0.1, 0.1, 1.229, 0.522, 0.1, 0.1,  0.517,...
%         0.504, 0.1, 0.101, 0.156, 0.547, 0.442, 0.59];
% %  Self Adaptive Harmony Search
% d_SAHS=[1.86, 0.521, 0.1, 0.1, 1.271, 0.509, 0.1, 0.1,...
%     0.485, 0.501, 0.1, 0.1, 0.168, 0.584, 0.433, 0.520];
% % Teaching-Learning-Based Optimization
% d_TLBO=[1.906, 0.506, 0.1, 0.1, 1.2617, 0.5111, 0.1, 0.1,...
%     0.5317, 0.5159, 0.1, 0.1, 0.156, 0.549, 0.409, 0.569];
%% a list of Reliability-based baseline designs (for the 72-truss)
% RBDOs in https://doi.org/10.1016/j.advengsoft.2015.11.001
% d_ICDE=[0.101, 0.533, 0.433, 0.488, 0.521, 0.508, 0.112, 0.105,  1.225, 0.521, 0.116, 0.111, 1.843, 0.491, 0.114, 0.118];
% d_NDL_ICDE=[0.1, 0.684, 0.525, 0.674, 0.720, 0.679, 0.10, 0.103,  1.631, 0.665, 0.10, 0.101, 2.471, 0.662, 0.1, 0.1];
% d_SORA_ICDE=[0.1, 0.682, 0.533, 0.661, 0.727, 0.670, 0.10, 0.10,  1.665, 0.668, 0.10, 0.10, 2.447, 0.667, 0.1, 0.1];
% T. Vo-Duy et al https://www.worldscientific.com/doi/abs/10.1142/S0219876219500166
 dn_Duy_RBMOO_E=[20.64, 9.59, 2.13, 0.65, 20.65, 10.32, ...
     0.77, 0.88, 11.02, 10.52, 0.65, 0.65, 2.39, 10.36, 3.45, 5.18]*0.155; % RBMOO - E % areas from  [cm^2]-->[in^2]
 % dn_Duy_RBMOO_C=[20.51, 9.57, 2.13, 0.65, 16.86, 10.2,....
%      0.65, 0.65, 10.44, 10.65, 0.65, 0.65, 4.11, 10.2, 0.65, 0.75]*0.155; % RBMOO - C % areas from  [cm^2]-->[in^2]
% % Truong et al .........https://doi.org/10.1016/j.advengsoft.2018.03.006
%  dn_RBDO_Truong_aeDE= [7419.3, 1858.06, 1045.1 1045.1 4658.05 1283.8 ...
%      1045.1 1045.15 2290.31 1045.1 1045.1 1045.1 1045.1 1045.1 1045.1 1045.1]*0.00155; % areas from  [mm^2]-->[in^2] (discrete case)
%  %% a list of Hybrid fuzzy/imprecise designs (for the 72-truss)
% d5_fuzzy=[0.13018 0.452 0.3297 0.4709 0.43942 0.4348  0.1 0.1 1.09388 0.43269 0.1 0.1 1.5946 0.43349 0.1, 0.1 ];%https://doi.org/10.1016/j.compstruc.2006.08.017
%
%% START!!!!!
load('Designs_72trussProblem.mat')
%1) compute nomilan desing reliability
for d=1:size(Designs_72trussProblem,2)
    dnominal_toCheck=Designs_72trussProblem(d).Areas;
    Rel_nominal{d}=ComputeReliabilityPerformance(dnominal_toCheck, delta,g_fun);
    Weight_nominal(d)=J_fun(dnominal_toCheck);
end
%% 2) run scenario program optimizer
VaR=0;
RHO=10^5; 
dn=dn_Duy_RBMOO_E;
tic
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',25000,'display','iter','UseParallel',true);
[Results,X_opt,J_opt,~,optresults]=ScenarioOptimizerCVaR_program1(VaR,w_fun,delta,J_fun,RHO,dn,LBd,UBd,options);
CompTime=toc;
Rel_sp_opt=ComputeReliabilityPerformance(Results.dopt ,delta_testing,g_fun);  % get reliability
sN=Results.Support.Size; % get size of the support sceanarios
 Weight_sp_opt=J_fun(Results.dopt);
 
% Robustness
beta=10^-4;
OutWeJ = getWaitandJudgeEpsilon(sN ,N,beta);% get Wait-and-judge reliability
epsilon_wait_and_judge =OutWeJ(end); % get Wait-and-judge reliability for the computed support size
[EpsilonLU(1), EpsilonLU(2)]= epsLU(sN,N, beta);
Epsilon_ind=Get_Epsilon_Individual_Requirements(g_fun,delta,Results.dopt,beta,VaR);

