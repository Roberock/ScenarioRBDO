function [g_fun,delta,dn,LBd,UBd]=Select_Case_Study(N,Case_LSF)
%% This function load the case study
% INPUT
% load N scenarios (interger N<=1e4)
% Case_LSF case study index 1,2,3,4,5,6
%% OUTPUT
% g_fun=performance function
% delta = the scenarios
% dn= nominal desing
% LBd lower bound on the desing values
% LBd upper bound on the desing values
addpath([pwd '\G_functions'])
addpath([pwd '\Scenarios'])
if Case_LSF==1 % case 1: 2 Linear LSFs with competitive failure domains
    %% Artifically created by Roberto Rocchetta 2019
    %% load scenarios
    TempLoad=load('Delta_N1e4_Grooteman2011Modified.mat');%N=1e4; % number of scenarios
    delta=TempLoad.delta(1:N,:);%  X1=normrnd(0,1,[N,1]);  X2=normrnd(0,2,[N,1]);  delta=[X1 X2]; % set of scenario samples
    %% load system parameters, g function
    g_fun=@(delta,d) g_LinearCompetitive(delta,d);
    a=1; b=1; c=1;
    dn=[a,b,c];
    LBd=[.5 .5 .5  ];UBd=[2 2 2 ];
elseif Case_LSF==2 % case 2: 2 non-Linear non-convex polynomial LSFs with competitive failure domains
    %% G function: LSF modified from [Grooteman2011]
    %% load scenarios
    TempLoad=load('Delta_N1e4_Grooteman2011Modified.mat');%N=1e4; % number of scenarios
    delta=TempLoad.delta(1:N,:);%  X1=normrnd(0,1,[N,1]);  X2=normrnd(0,2,[N,1]);  delta=[X1 X2]; % set of scenario samples
    %% load system parameters, g function
    dn=[2.5,0.2 ,0.06]; % nominal desing vector
    g_fun=@(delta,d) g_Grooteman2011_modified(delta,d);
    LBd=[0.5, -2, -0.3 ];% Assume lower bound on the desing variables and gamma
    UBd=[4,    2,  0.3 ];% Assume lower bound on the desing variables and gamma
elseif Case_LSF==3 % to be defined (e.g. high Ng number of narbitrary (e.g nom polynomial) and competitive LSFs)
    %%  Parallel System Problem 4 g
    %% load scenarios
    TempLoad=load('Delta_N1e4_Grooteman2011Modified.mat');%N=1e4; % number of scenarios
    delta=TempLoad.delta(1:N,:);%  X1=normrnd(0,1,[N,1]);  X2=normrnd(0,2,[N,1]);  delta=[X1 X2]; % set of scenario samples
    %% load system parameters, g function
    dn=[0.1,2,1,3,3]; % nominal desing vector
    % dn=[0.1,2,1,3.5,3]; % nominal desing vector
    g_fun=@(delta,d) g_4D_parallelSystem(delta,d);
    LBd=[0.05, 1.5, 0.5,2,2];% Assume lower bound on the desing variables and gamma
    UBd=[0.15, 2.5, 1.5,4,4];% Assume lower bound on the desing variables and gamma
elseif Case_LSF==5 % to be defined (e.g. high Ng number of narbitrary (e.g nom polynomial) and competitive LSFs)
    %%  Antenna Tower Case Study (25 beams bending tower)
    %% load scenarios
    % E=normrnd(1e7,0.05*1e7,[N,25]);
    % theta=unifrnd(-5/180*pi,5/180*pi,[N,1]);
    % phi=unifrnd(-pi,+pi,[N,1]);
    %  delta=[E,theta,phi];
    TempLoad=load('delta_Antenna_N1e4''delta.mat');
    delta=TempLoad.delta(1:N,:);
    %% load system parameters, g function
    d1=0.4;  d2=0.1;
    d3=3.4;  d4=1.3;
    d5=0.9;  d6=1.0;
    dn=[d1 d2 d3 d4 d5 d6]; % nominal desing vector
    g_fun=@(delta,d) g_AntennaTower(delta,d);
    LBd=0.8.*[dn];% Assume lower bound on the desing variables and gamma
    UBd=1.2.*[dn];% Assume lower bound on the desing variables and gamma
    %     BeamLengths=[75.0000  130.5038  130.5038  130.5038  130.5038  106.8000  106.8000  106.8000  106.8000   75.0000   75.0000   75.0000   75.0000  181.1422 181.1422  181.1422  181.1422  181.1422  181.1422  181.1422  181.1422  133.4635  133.4635  133.4635  133.4635];
    %     fun_volume=@(d,BeamLengths)  beamLengths' .* d;
elseif Case_LSF==4 % SPEED REDUCER BOX 6G 5 random variables 6 design variables
    %% load scenarios 
    %  mu=[4.7835e+03, 1.3977e+06, 1.93, 16.9 * 10^6 ,157.5*1e6];
    % Cov=[0.15, 0.1, 0.15, 0.15, 0.15];
%     delta(:,1)=normrnd(4.7835e+03,4.7835e+03.*0.15,[N,1]);
%     delta(:,2)=normrnd(1.3977e+06,1.3977e+06.*0.1,[N,1]);
%     delta(:,3)=normrnd(1.93,1.93*0.15,[N,1]);
%     delta(:,4)=normrnd(16.9 * 10^6,16.9 * 10^6*0.15,[N,1]);
%     delta(:,5)=normrnd(157.5* 10^6,157.5* 10^6*0.15,[N,1]) ;
%     delta(:,6)=unifrnd(5.1,5.2,[N,1]) ;
    TempLoad=load('delta_SpeedReducer.mat');
    delta=TempLoad.delta(1:N,:);
    %% load system parameters, g function
    dn=[3, 0.75, 7.5, 7.9, 3 ];% d6 ]; % nominal desing vector
    g_fun=@(delta,d) g_SpeedReducer(delta,d);
    LBd=[2.6, 0.7, 7.3, 7.3, 2.9 ];% Assume lower bound on the desing variables and gamma
    UBd=[3.6, 0.8, 8.3 ,8.3, 3.9 ];% Assume lower bound on the desing variables and gamma
elseif Case_LSF==6 % Crespo problem (benchmark robust control)
    %% load scenarios
    addpath([pwd '\G_functions\ContollerCrespo'])
    TempLoad=load('ScenarioOptData_UncertaintyLargeonM1M2.mat');
    delta=TempLoad.ScenarioOptData.delta(1:N,:); %N<=1e3;
    %% load system parameters, g function
    dn=[-0.1324 0.3533 0.6005 0.0728 0.5503 1.4175 2.6531 2.4802 1];
    g_fun=@(delta,d) g_ControllerCrespo(delta,d);
    LBd=dn-abs(dn).*0.8;% Assume lower bound on the desing variables and gamma
    UBd=dn+abs(dn).*0.8;% Assume lower bound on the desing variables and gamma
    
elseif Case_LSF==7     % modified Rosenbrock FUNCTION 50 random variables 10 desing variables
    
    prompt = 'how many desing variables? ';
    Nd = input(prompt);
    %% load scenarios
    TempLoad=load('ScenarioOptData_Rosenbrock2.mat');
    delta=TempLoad.delta(1:N,:); %N<=1e4;
    %% load system parameters, g function
    dn=3.*ones(1,Nd); % initial desing in 0;
    g_fun=@(delta,d) g_Rosenbrock(delta,d);
    LBd=-30.*ones(1,Nd) ;% Assume lower bound on the desing variables
    UBd=+30.*ones(1,Nd) ;% Assume lower bound on the desing variables
    
 elseif Case_LSF==8     % modified Rosenbrock FUNCTION 50 random variables 10 desing variables
 
    %% load scenarios
    TempLoad=load('Delta_N1e4_Grooteman2011Modified.mat');
    delta=TempLoad.delta(1:N,:); %N<=1e4;
    %% load system parameters, g function
    Nd=2;
    dn=ones(1,2); % initial desing in 1;
    g_fun1=@(delta,d) +d(1).*sin(d(1)*(delta(:,2)));
    g_fun2=@(delta,d) -d(1).^5*(abs(delta(:,1)))+1/d(2)*(delta(:,2)) ;
    g_fun=@(delta,d) [g_fun1(delta,d) g_fun2(delta,d)];
    LBd=-5.*ones(1,Nd) ;% Assume lower bound on the desing variables
    UBd=+5.*ones(1,Nd) ;% Assume lower bound on the desing variables
       
end

end