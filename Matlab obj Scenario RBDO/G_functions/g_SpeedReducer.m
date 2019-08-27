function g=g_SpeedReducer(delta,d)
% MODIFIED Speed Reducer Design Problem
% see e.g. Design Optimization of a Speed Reducer Using Deterministic Techniques Ming-Hua Lin,1 Jung-Fa Tsai,2 Nian-Ze Hu,3 and Shu-Chuan Chang2,4
% x1: the face width
% x2: the module of the teeth
% x3: the number of teeth on pinion
% x4: the length of the first shaft between bearings 
% x5: the length of the second shaft between bearings
% x6: diameter of the first shaft
% x7: and the diameter of the second shaft
%% design variables
x1=d(1); 
x2=d(2);%  
x3=17; % assum integer fixed to 17
x4=d(3);% 
x5=d(4);
x6=d(5); 
%x7=d(6);
%% delta Lets assume we have already a set of shaft manufactured in a certain way
% uncertainty is affecting x1,x4,x5 due to variability in the manufacturing process
% delta(:,1)=2*M*q/k_g (originally==27)   where M=Bending moment?    q=tooth form factor 2.54   k_g=Permissible bending stress of gear teeth 900
% M=4.7835e+03
% delta(:,2)=2*B*M/P_d^2 (originally=397.5) where M=Bending moment?   B= elastic coefficient dependent on the modulus of elasticity  P_d=Permissible surface compressive stress 5800 kG cm -2
% B=1.3977e+06 
% delta(:,3)= related to the force P transmitted (originally=1.93)  

% delta(:,4)= (originally=16.9 × 10^6)  relted to ...Mz1/Wx1    <= 1100 kGcm-2... Mzi bending and torsional moments for shafts i 

% delta(:,5)= (originally=157.5*1e6)  relted to ...Mz2/Wx2    <= 850 kGcm-2...   Wx1 =strength section modulus
% delta(:,6)=diameter of the second shaft 

%% G1
kg=900;q=2.54;
M=delta(:,1); % Bending moment
g1=(2*M.*q./kg)./(x1.*x2.^2.*x3)-1; % Bending condition i.e. actual bending stress of gear teeth 
%% G2
P_d=5800;
B=delta(:,2);% B= elastic coefficient dependent on the modulus of elasticity
g2=(2.*B.*M./P_d.^2)./(x1.*x2.^2.*x3.^2)-1; % compressing stress limitation 
%% G3-9
g3=delta(:,3).*x4.^3./(x1.*x3.*x6^4)-1;%  Transverse deflection of shaft 1
g4=delta(:,3).*x5.^3./(x2.*x3.*delta(:,6).^4)-1;% Transverse deflection of shaft 2 
g5=1./(110.*x6^3).*sqrt(((745.*x4)./(x2.*x3)).^2+delta(:,4))-1; %Substitute stress condition for shaft 1    Mz1/Wx1    <= 1100 kG cm-2" 
g6=1./(85.*delta(:,6).^3).*sqrt(((745.*x5)./(x2.*x3)).^2+delta(:,5))-1; % Substitute stress condition for shaft 2   Mz2/Wx2    <= 850 kG cm-2" 
% g7=x2.*x3./40-1; % Overall dimensions condition x2(x3+x3)< 160. 
% g8=5*x2./x1-1;  %Relative face  width conditions 1 x1/x2>=5
% g9=x1./(12.*x2)-1; %Relative face  width conditions 2 x1/x2<=12
% g10=(1.5*x6+1.9)./x4-1; % Design condition 1..... 1.5x6+ 1.9 <=x4  
% g11=(1.1*x7+1.9)./x5-1; % Design condition 2 ....1.1x7+ 1.9 <=x5  
%%
g=-[-g1 -g2 -g3 -g4 -g5 -g6];% 

end
