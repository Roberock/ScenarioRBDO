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
    LBd=[0.05, 1.5, 0.5, 2, 2];% Assume lower bound on the desing variables and gamma
    UBd=[0.15, 2.5, 1.5, 4, 4];% Assume lower bound on the desing variables and gamma

elseif Case_LSF==5    % Roberto s Case STudy non-Convex
    %% generate scenarios
    rng default
    delta(:,1)=normrnd(0.5,1.5,[N,1]);
    delta(:,2)=sqrt(unifrnd(1,4,[N,1])+delta(:,1).^2 +delta(:,1));
    %% nominal desing
    dn= [1 0.5];
    %% G1:
    g_fun1=@(delta,d)   +0.5+cos(d(1)^2+d(1)*d(2))*(delta(:,2)-delta(:,1));%LSF linear dependency with d
    %% G2:
    g_fun2=@(delta,d)  (-d(2)*d(1).*delta(:,1).^2).*sin(d(2)+delta(:,2)-delta(:,1))./(1+sin(d(2)*d(1)+delta(:,2).^2-delta(:,1)).^2) ;
    g_fun=@(delta,d) [g_fun1(delta,d) g_fun2(delta,d) ]; %g>0 is failure
    %% lower upper bounds
    LBd=[1  -1 ];% Assume lower bound on the desing variables
    UBd=[3 1];% Assume lower bound on the desing variables
    elseif Case_LSF==4 % SPEED REDUCER BOX 6G 5 random variables 6 design variables
    %% load scenarios
    %  mu=[4.7835e+03, 1.3977e+06, 1.93, 16.9 * 10^6 ,157.5*1e6];
    % Cov=[0.15, 0.1, 0.15, 0.15, 0.15];  delta(:,2)=normrnd(1.3977e+06,1.3977e+06.*0.1,[N,1]);
    %     delta(:,3)=normrnd(1.93,1.93*0.15,[N,1]);   delta(:,4)=normrnd(16.9 * 10^6,16.9 * 10^6*0.15,[N,1]);
    %     delta(:,5)=normrnd(157.5* 10^6,157.5* 10^6*0.15,[N,1]) ;   delta(:,6)=unifrnd(5.1,5.2,[N,1]) ;
    TempLoad=load('delta_SpeedReducer.mat');  delta=TempLoad.delta(1:N,:);
    %% load system parameters, g function
    dn=[3, 0.75, 7.5, 7.9, 3 ];% d6 ]; % nominal desing vector
    g_fun=@(delta,d) g_SpeedReducer(delta,d);
    LBd=[2.6, 0.7, 7.3, 7.3, 2.9 ];% Assume lower bound on the desing variables and gamma
    UBd=[3.6, 0.8, 8.3 ,8.3, 3.9 ];% Assume lower bound on the desing variables and gamma
elseif Case_LSF==6 %  Canti Leaver Beam with non-random length
    % http://web.mst.edu/~dux/repository/me360/Examples/Examples_Part%203/3_1.pdf
    %%   % minimize (d(1)*d(2)*d(3))*(1+gamma)
    %  P=(Px,Py,S,E)
    rng default
    delta =[normrnd(500,100^2,[N,1]),...
        normrnd(500,100^2,[N,1]),...
        normrnd(40*10^3,2*10^6,[N,1]),...
        normrnd(29*10^6,3*10^6,[N,1]),...
        normrnd(0,0.1^2,[N,1])] ;
    
    %% G function: cantileaver
    g_fun=@(delta,d)  G_cantileaver_fun(delta,d) ;
    
    %% load system parameters, g function
    dn = [2.0,3.0,40.0]; % dvs=[b,h,Length] nominal
    LBd = [1.0,1.0,50.0]; % lower bounds of design
    UBd = [20.0,20.0,100.0]; % upper bounds of dvs
     
elseif Case_LSF==7 % Crespo problem (benchmark robust control)
    %% load scenarios
    addpath([pwd '\G_functions\ContollerCrespo'])
    TempLoad=load('ScenarioOptData_UncertaintyLargeonM1M2.mat');
    delta=TempLoad.ScenarioOptData.delta(1:N,:); %N<=1e3;
    %% load system parameters, g function
    dn=[-0.1324 0.3533 0.6005 0.0728 0.5503 1.4175 2.6531 2.4802 1];
    g_fun=@(delta,d) g_ControllerCrespo(delta,d);
    LBd=dn-abs(dn).*0.8;% Assume lower bound on the desing variables and gamma
    UBd=dn+abs(dn).*0.8;% Assume lower bound on the desing variables and gamma
    
elseif Case_LSF==8     % modified Rosenbrock FUNCTION 50 random variables 10 desing variables
    
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
 elseif Case_LSF==9 % to be defined (e.g. high Ng number of narbitrary (e.g nom polynomial) and competitive LSFs)
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
       
% elseif Case_LSF==10    %
%     
%     %% load scenarios
%     TempLoad=load('delta_D31cluster.mat');
%     delta=TempLoad.delta_D31cluster(1:N,:);
%     delta(:,1)=delta(:,1);delta(:,2)=delta(:,2);
%     %% load system parameters, g function
%     Nd=6; dn= [1,1,1,1,1,1]; % initial desing in 1; 
%     %% G1:
%     [g_fun] =@(delta,d) Gfun_hart3(delta,d);
%     UBd=1.*ones(1,Nd) ;% Assume lower bound on the desing variables
%     LBd(end)=0.*ones(1,Nd);
%     
    
% elseif  Case_LSF==10
%     %% load scenarios
%     TempLoad=load('delta_D31cluster.mat');
%     delta=TempLoad.delta_D31cluster(1:N,:);
%     delta(:,1)=delta(:,1);delta(:,2)=delta(:,2);
%     %% load system parameters, g function
%     Nd=6;
%     dn= [1,1,1,1,1,1]; % the desing variables are he
%     
%     Emodulus=188e9; % Modulus of elasticity over Lentgh(207e9/1.1)
%     
%     kmatrix=[1.454 0.454 -1 0 0 0 0 0 -0.454 -0.454 0 0 ;% K is the total stiffness matrix
%         0.454 0.454 0 0 0 0 0 0 -0.454 -0.454 0 0;
%         -1 0 2.908 0 -1 0 -0.454 -0.454 0 0 -0.454 0.454;
%         0 0 0 1.908 0 0 -0.454 -0.454 0 -1 0.454 -0.454;
%         0 0 -1 0 1.454 -0.454 0 0 -0.454 0.454 0 0;
%         0 0 0 0 -0.454 1.454 0 -1 0.454 -0.454 0 0;
%         0 0 -0.454 -0.454 0 0 1.454 0.454 -1 0 0 -1;
%         0 0 -0.454 -0.454 0 -1 0.454 1.454 0 0 0 0;
%         -0.454 -0.454 0 0 -0.454 0.454 -1 0 2.908 0 -1 0;
%         -0.454 -0.454 0 -1 0.454 -0.454 0 0 0 1.908 0 0;
%         0 0 -0.454 0.454 0 0 0 0 -1 0 1.454 -0.454;
%         0 0 0.454 -0.454 0 0 0 0 0 0 -0.454 0.454];
%     
%     
%     DisPlacement=[0 0 d3 d4 d5 d6 d7 d8 d9 d10 0 0]';
%     Force=[f1 f2 f3 d4 0 -5000 f7 f8 f9 f10 f11 f12]';
%     
%     Totalmatrix=@(x) Emodulus*d *kmatrix*DisPlacement;
%     
%     c(1)=DisPlacement-0.05; % maximum allowed displacement is less than 5 cm
%     ceq(1)=Totalmatrix-Force;
%     
end

end