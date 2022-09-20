function [g_fun,delta,dn,LBd,UBd,DGM,Nd,Ng]=Select_Convex_Case_Study(N,Case_LSF)
%% This function load the case study
%% INPUT
% N= number of scenarios (interger >0)
% Case_LSF= the case study, integer 1,2,3...
%% OUTPUT
% g_fun = performance functions
% delta = the scenarios with the uncertain parameters
% dn = baseline desing
% Nd = number of desing variables
% LBd = lower bound on the desing values
% LBd = upper bound on the desing values
% DGM = the data-generating-mechanism (for testing)
%% List of cases
% case 1: Two requirements non-linear in delta (modified from [Grooteman2011])
% case 2: Four requirements non-linear in delta (modified from [Shueller2004])
% case 3: Find engingeering or portfolio optimization example.......

%% START
addpath([pwd '\G_functions'])
rng default % rest random seed for reproducibility

if Case_LSF==0 %   linear competitive
    
    %% Define Data-Generating-Mechanism and sample N scenarios delta
    DGM= @(N)  mvnrnd([0 0],[1  0; 0 2^2],N); % data generating mechanism
    delta=DGM(N); % available samples
    %% Define Data-Generating-Mechanism and sample N scenarios delta
    dn=[1,1,1];
    Nd=length(dn); Ng=2;
    LBd= [0.5, 0.5, 0.5] ;  % design's lower bounds
    UBd= [2 , 2 , 2];  % design's upper bounds
    g1=@(d,delta) delta(:,2)./d(1) + delta(:,1)./d(2) - d(3);
    g2=@(d,delta) d(1)*(delta(:,1)) - delta(:,2)./d(2) - d(3);
    g_fun=@(d,delta) [g1(d,delta) g2(d,delta)];
    
elseif Case_LSF==1 %    Grooteman 2011
    %% Define Data-Generating-Mechanism and sample N scenarios delta
    rng default % rest random seed for reproducibility
    DGM= @(N)  mvnrnd([0 0],[1  0; 0 2^2],N); % data generating mechanism
    delta=DGM(N); % available samples
    %% Define Data-Generating-Mechanism and sample N scenarios delta
    dn=[2.5,0.2,0.06];
    Nd=length(dn);
    LBd= [0.5, -2, -0.3] ;  % design's lower bounds
    UBd= [4 , 2 ,  0.3];  % design's upper bounds
    g_fun=@(d,delta) g_Grooteman2011_modified(delta,d);
    Ng=length(g_fun(dn,delta(1,:)));
elseif Case_LSF==2 % from Shueller 2004
    %% Define Data-Generating-Mechanism and sample N scenarios delta
    rng default % rest random seed for reproducibility
    DGM= @(N)  mvnrnd([0 0],[1.2  -0.9; -0.9 1.2],N); % data generating mechanism
    delta=DGM(N); % available samples
    
    %% g function, baseline desing and desing bounds
    dn=[0.2,0.8801,1,6];
    Nd=length(dn);
    LBd= [-0.5, 0.1, 1, 5] ;  % design's lower bounds
    UBd= [0.5 , 2 ,  2, 7];  % design's upper bounds
    g_fun=@(d,delta) G_fun_Shueller2004(d,delta);
    Ng=length(g_fun(dn,delta(1,:)));
    
elseif Case_LSF==3  %
    %[] B. J. Bichon et al. , Efficient Global Reliability Analysis for Nonlinear Implicit Performance Functions, AIAA JOURNAL 2008
    %% Define Data-Generating-Mechanism and sample N scenarios delta
    rng default % rest random seed for reproducibility
    DGM= @(N)  mvnrnd([1.5 2.5],[1  0;0 1],N); % data generating mechanism
    delta=DGM(N); % available samples
    %% g function, baseline desing and desing bounds
    dn=[2,1];
    Nd=length(dn);
    LBd= [1, -1.5 ] ;  % design's lower bounds
    UBd= [3 , 1.5  ];  % design's upper bounds
    g_fun=@(d,delta) ((delta(:,1).^2+d(1)^2).*(delta(:,2)-d(2)))./20-d(2).*sin(5.*delta(:,1)./2)-2;
    Ng=length(g_fun(dn,delta(1,:)));
elseif Case_LSF==4  % manufacturer produces goods problem
    % Authors:  Garatti Campi 2019, Riskandcomplexityinscenariooptimization, MathematicalProgramming https://doi.org/10.1007/s10107-019-01446-4
    % modified by R. Rocchetta 2020
    %% Define Data-Generating-Mechanism and sample N scenarios delta
    rng default % rest random seed for reproducibility
    Nd=50; % Nd workplaces
    Nk=5;  % Number of resources k used to produce the
    alphas= @(N) unifrnd(10,50,[N,Nk]);
    gammas= @(N)  unifrnd(-2.5,+2.5,[N,Nk*Nd]);
    DGM= @(N) [alphas(N) gammas(N)];
    delta =DGM(N);
    %% g function, baseline desing and desing bounds
    dn=0.01*ones(1,Nd); % quantity of resource to be produced in the Nd workplaces
    LBd= [zeros(1,Nd)] ;  % design's lower bounds
    UBd= [10*ones(1,Nd)];    % design's upper bounds
    a=ones(1,Nk); % treshold on the workplace budget (?)
    a=[1 3 1 5 4];
    g_fun= @(d,delta) g_manufacturer_goods(d,delta,Nk,a);
    Ng=Nk;
elseif Case_LSF==5 %     % OTL CIRCUIT FUNCTION (non-convex..)
    % Authors: Sonja Surjanovic, Simon Fraser University
    %          Derek Bingham, Simon Fraser University
    % modified by Rocchetta 2020
    %% Define Data-Generating-Mechanism and sample N scenarios delta 
    DGM= @(N) mvnrnd([0,0,0,0,0,0],[7.5,3.5,0.15,0.125,0.06,15].^2,N);
    delta =DGM(N);
    %% g function, baseline desing and desing bounds
    dn=[ 100 , 47.5  , 1.75  , 1.85 ,  0.7250 , 175];
    Nd=length(dn);
    LBd= [50,25,0.5,1.2,0.25,50] ;  % design's lower bounds
    UBd= [150,70,3,2.5,1.2,300];  % design's upper bounds
    g_fun= @(d,delta) gfun_otlcircuit(d,delta);
    Ng=length(g_fun(dn,delta(1,:)));
    
    
elseif Case_LSF==15 % big size delta
    %% Define Data-Generating-Mechanism and sample N scenarios delta 
    DGM= @(N) betarnd(2,0.2,[1,N])+ betarnd(0.1,2,[1,N]).^2;
    delta =DGM(N);
    %% g function, baseline desing and desing bounds 
    yfun = @(delta,theta) theta(3)*delta(:,1).^2 .*sqrt(4*delta(:,2)+0.25*delta(:,3)+1)+(sum(theta(1:2))-0.5).^2; % system function
    g_fun= @(d,delta) yfun(delta,d)-3; %  
    dn=[1 2 3];
    Nd=length(dn);
    LBd= [-10 -10 -10] ;  % design's lower bounds
    UBd= [10 10 10];  % design's upper bounds  
    
elseif Case_LSF==15 % big size delta
    Ndelta=90;
    Nd=6;
    %% Define Data-Generating-Mechanism and sample N scenarios delta 
    Lb_delta = ones(1,Ndelta);
    Ub_delta = ones(1,Ndelta).*(randi(10,[1,Ndelta])+1);
    DGM=@(N) unifrnd(Lb_delta,Ub_delta,N);
    delta=zeros(N,Ndelta);
    for k=1:N
        delta(k,:)=unifrnd(Lb_delta,Ub_delta); % available samples
    end
    %% g function, baseline desing and desing bounds
    dn=ones(1,Nd);
    LBd= dn-0.5;  % design's lower bounds
    UBd= dn+0.5;  % design's upper bounds
    % G1
    g_1=@(d,delta)  66 - prod(delta(:,81:83),2)*d(1).*(sum(sin(delta(:,1:15))+d(5).*delta(:,16:30),2));
    %G2
    g_2=@(d,delta)  72 -prod(delta(:,84:86),2)*d(4)*d(2).*(sum(sin(delta(:,31:45))+d(5).*delta(:,46:60),2));
    %G3
    g_3=@(d,delta)  58- prod(delta(:,87:90),2)*d(3).*(sum(sin(delta(:,61:70))-d(1)*d(5).*delta(:,71:80),2));
    %% G all
    g_fun=@(d,delta) [g_1(d,delta)  g_2(d,delta)  g_3(d,delta) ];
    Ng=length(g_fun(dn,delta(1,:)));
    %% ALTERNATIVE CASE STUDIES
    
    %% %%%%%%%%%%%%%%%%%%%   problem 1: Linear-Competitive (g1,g2) %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ng=2; % number 2of performance functions
    % Nd=2; % number of desing variables
    %g1= @(d,delta)  -delta(:,1).*d(1) + d(2) +delta(:,2) ;   % linear Performance Functions 1
    %g2= @(d,delta)  -delta(:,2).*d(2) + d(1) +delta(:,1);  % linear Performance Functions 2
    %g_fun=@(d,delta)[g1(d,delta), g2(d,delta) ]; % Performance Functions
    % dn=[2,3]; % nominal desing
    % LBd=[-5,-5]; UBd=[0,0];  % design bounds
    % J_fun= @(d)  sum(d); % first part of the cost function J(d)
    %% %%%%%%%%%%%%%%%%%%%  problem 3: Roof Trust reliability problem %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  delta(:,1)=normrnd(20000,1400,[N,1]);   %  mu=20000 ; std=1400; %normal random Load
    %  delta(:,2)=normrnd(12,0.12,[N,1]);    %  mu=12 ; std=0.12; %normal random Length factors for truss members
    %  delta(:,3)=normrnd(10^11,6*10^9,[N,1]);   %  mu=10^11 ; std=6*10^9; %Elastic modulus of steel
    %  delta(:,4)=normrnd(2*10^10,1.2*10^9,[N,1]);   %  mu=2*10^10 ; std=1.2*10^9; %Elastic modulus of concrete
    %  delta(:,5)=normrnd(0,5.9853*10^(-5),[N,1]); %  mu=0 ; std=5.9853*10^(-5); %error on the cross-sectional area of steel bars
    %  delta(:,6)=normrnd(0, 0.0048,[N,1]); %  mu=0 ; std= 0.0048; %error on the cross-sectional area of concrete bars
    % g_fun=@(d,delta) G_roofTrust(d,delta);
    % Ng=1; % number of performance functions
    % Nd=2; % number of desing variables
    % dn=[0.001, 0.042]; % nominal desing
    % LBd=[ 0.0006,0.018 ];
    % UBd=[0.0012,0.063];
    % J_fun=@(d) 20224*d(1)+364*d(2);
    % Reference Solutions at: https://arxiv.org/pdf/1904.11424.pdf THRESHOLD SHIFT METHOD FOR RELIABILITY-BASED DESIGN
    %opt-1: dopt1=[ 10.1*10^-4, 3.45*10^-2]; Pf(d*)= ; J(d*)=32.97;% BetaPf=1.989; % fun call 1702
    %opt-2: dopt2=[6*10^-4, 3*10^-2]; Pf(d*)=; J(d*)= 23.0; BetaPf= 1.65;  % fun call 5,224
    %% %%% %%%%%  problem 4 :   %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ng=12;Nd=11;
    % delta=[normrnd(70*10^9,10^7,[N,11]),normrnd(-10^6,2*10^4,[N,4])]; delta(:,14)=-delta(:,14);
    % g_fun=@(d,delta) g_2Dtrussproblem(d,delta); % 11 elements truss problem
    % dn=[ones(1,11)*3*10^(-4)];
    % LBd=[ones(1,11)*2*10^(-4)];
    % UBd=[ones(1,11)*4*10^(-4)];
    % L=[3,3,4,3,4,3,4,4,3,3,3]; % length of the elements
    % J =@(d) sum(L.*d); % cost proportional to the volume
    %% %%%%%%%%%%%%%%%%%%%  problem 5 :   %%%%%%%%%%%%%%%%%%%%%%%%%%
    % addpath('C:\Users\Roberto.Rocchetta\Desktop\NASA-NIA\ESREL 2020\scenario theory\G_functions')
    % Ng=5; % number of performance functions
    % Nd=24; % number of desing variables
    % Ndelta=12;
    % delta=normrnd(0,1,[N,Ndelta]);
    % g_fun=@(d,delta) G_function_controller(d,delta);
    %  Ng=size(g,2);
    % dn=[ones(1,Nd)]; % nominal desing
    % LBd=[-2*ones(1,Nd)]; UBd=[2*ones(1,Nd)];  % design bounds
    % J_fun= @(d) trace(reshape(d(1:16),4,4)); % first part of the cost function J(d)
    
    
    
end

end